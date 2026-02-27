# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import os
import sys
from pathlib import Path
import numpy as np
from importlib import resources
from ctypes import (
    CDLL,
    CFUNCTYPE,
    POINTER,
    byref,
    string_at,
    cast,
    addressof,
    memmove,
    c_double,
    c_int64,
    c_int32,
    c_bool,
    c_void_p,
    c_char_p,
    c_char,
    Structure,
)
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any


# load the opentrustregion library, fallback to testsuite in case opentrustregion was
# statically compiled
if sys.platform == "darwin":
    ext = "dylib"
elif sys.platform == "win32":
    ext = "dll"
else:
    ext = "so"
lib_candidates = [
    f"libopentrustregion.{ext}",
    f"libopentrustregion_32.{ext}",
    f"libopentrustregion_64.{ext}",
]
if sys.platform != "win32":
    # libopentrustregion symbols not exported through libotrtestsuite on windows
    lib_candidates.append(f"libotrtestsuite.{ext}")
lib = None

# try to load from installed package (conda)
conda_prefix = Path("/opt/anaconda1anaconda2anaconda3")
for lib_name in lib_candidates:
    conda_path_unix = conda_prefix / "lib" / lib_name
    conda_path_wind = conda_prefix / "Library" / "bin" / lib_name
    if conda_path_unix.exists():
        lib = CDLL(str(conda_path_unix))
    elif conda_path_wind.exists():
        lib = CDLL(str(conda_path_wind))

# fallback: try to load from installed package (site-packages)
if lib is None:
    for lib_name in lib_candidates:
        lib_path = resources.files("pyopentrustregion") / lib_name
        if lib_path.is_file():
            lib = CDLL(str(lib_path))

# fallback: try to load from same directory (editable install)
if lib is None:
    for lib_name in lib_candidates:
        local_path = os.path.join(os.path.dirname(__file__), lib_name)
        if os.path.exists(local_path):
            lib = CDLL(local_path)

# fallback: try to load from ../build (development build)
if lib is None:
    for lib_name in lib_candidates:
        build_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../build", lib_name)
        )
        if os.path.exists(build_path):
            lib = CDLL(build_path)

# if all failed
if lib is None:
    raise FileNotFoundError(
        f"Cannot find either opentrustregion or testsuite library ({lib_candidates})"
    )

# determine integer size used in library
ilp64 = c_bool.in_dll(lib, "ilp64")
if ilp64.value:
    c_int = c_int64
else:
    c_int = c_int32

# define real type
c_real = c_double

# fixed size strings for keywords
kw_len = c_int.in_dll(lib, "kw_len_c").value

# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
hess_x_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
update_orbs_interface_type = CFUNCTYPE(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(hess_x_interface_type),
)
obj_func_interface_type = CFUNCTYPE(c_int, POINTER(c_real), POINTER(c_real))
precond_interface_type = CFUNCTYPE(
    c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)
project_interface_type = CFUNCTYPE(c_int, POINTER(c_real))
conv_check_interface_type = CFUNCTYPE(c_int, POINTER(c_bool))
logger_interface_type = CFUNCTYPE(None, c_char_p)


# define interface factories
def hess_x_interface_factory(
    hess_x: Callable[[np.ndarray, np.ndarray], None], n_param: int
) -> Any:
    """
    this function is a factory for the Hessian linear transformation interface
    """

    @hess_x_interface_type
    def hess_x_interface(x_ptr, hx_ptr):
        """
        this function interfaces the Hessian linear transformation
        """
        # convert trial vector pointer to numpy array
        x = np.ctypeslib.as_array(x_ptr, shape=(n_param,))
        hx = np.ctypeslib.as_array(hx_ptr, shape=(n_param,))

        # perform linear transformation
        try:
            hess_x(x, hx)
        except RuntimeError:
            return 1

        return 0

    # attach the Hessian-vector product function to the returned interface so that it
    # persists in Python to ensure that it is not garbage collected when the factory
    # completes
    hess_x_interface.hess_x = hess_x

    return hess_x_interface


def precond_interface_factory(
    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]], n_param: int
) -> Any:
    """
    this function is a factory for the preconditioning interface
    """
    if precond is None:
        return None

    @precond_interface_type
    def precond_interface(residual_ptr, mu_ptr, precond_residual_ptr):
        """
        this function interfaces the preconditioner
        """
        # convert pointers to numpy array and float
        residual = np.ctypeslib.as_array(residual_ptr, shape=(n_param,))
        mu = mu_ptr[0]
        precond_residual = np.ctypeslib.as_array(precond_residual_ptr, shape=(n_param,))

        # call preconditioner
        try:
            precond(residual, mu, precond_residual)
        except RuntimeError:
            return 1

        return 0

    return precond_interface


def project_interface_factory(
    project: Optional[Callable[[np.ndarray], None]], n_param: int
) -> Any:
    """
    this function is a factory for the projection interface
    """
    if project is None:
        return None

    @project_interface_type
    def project_interface(vector_ptr):
        """
        this function interfaces the projection function
        """
        # convert pointers to numpy array and float
        vector = np.ctypeslib.as_array(vector_ptr, shape=(n_param,))

        # call projection function
        try:
            project(vector)
        except RuntimeError:
            return 1

        return 0

    return project_interface


def conv_check_interface_factory(conv_check: Optional[Callable[[], bool]]) -> Any:
    """
    this function is a factory for the convergence check interface
    """
    if conv_check is None:
        return None

    @conv_check_interface_type
    def conv_check_interface(conv_ptr):
        """
        this function interfaces the convergence check
        """
        # call convergence check
        try:
            conv_ptr[0] = conv_check()
        except RuntimeError:
            return 1

        return 0

    return conv_check_interface


def logger_interface_factory(logger: Optional[Callable[[str], None]]) -> Any:
    """
    this function is a factory for the logging interface
    """
    if logger is None:
        return None

    @logger_interface_type
    def logger_interface(message):
        """
        this function interfaces the logging
        """
        # call logger
        logger(string_at(message).decode("utf-8"))

    return logger_interface


# define classes corresponding to C structs for settings
class SolverSettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("project", c_void_p),
        ("conv_check", c_void_p),
        ("logger", c_void_p),
        ("stability", c_bool),
        ("line_search", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_real),
        ("start_trust_radius", c_real),
        ("global_red_factor", c_real),
        ("local_red_factor", c_real),
        ("n_random_trial_vectors", c_int),
        ("n_macro", c_int),
        ("n_micro", c_int),
        ("jacobi_davidson_start", c_int),
        ("seed", c_int),
        ("verbose", c_int),
        ("subsystem_solver", c_char * (kw_len + 1)),
    ]


class StabilitySettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("project", c_void_p),
        ("logger", c_void_p),
        ("initialized", c_bool),
        ("conv_tol", c_real),
        ("n_random_trial_vectors", c_int),
        ("n_iter", c_int),
        ("jacobi_davidson_start", c_int),
        ("seed", c_int),
        ("verbose", c_int),
        ("diag_solver", c_char * (kw_len + 1)),
    ]


# define classes for settings
class Settings:

    c_struct: type[Structure]
    init_c_struct: Any

    def __init__(self):
        """
        this function initializes the settings class
        """
        # specify C interface signature
        self.init_c_struct.argtypes = [POINTER(self.c_struct)]
        self.init_c_struct.restype = None

        # call C-side initialization to populate defaults
        self.settings_c = self.c_struct()
        self.init_c_struct(byref(self.settings_c))

        # initializes all optional function pointers to None
        for field_info in self.settings_c._fields_:
            field_name, field_type = field_info[:2]
            if field_type is c_void_p:
                setattr(self, field_name, None)

    def set_optional_callback(self, attr_name, func_interface):
        """
        this function sets a callback in the C struct from a Python function while also
        keeping the interface alive in the Python object
        """
        # create c_void_p pointer, points to NULL if func_interface is None
        setattr(self.settings_c, attr_name, cast(func_interface, c_void_p))

        # keep interface alive
        setattr(self.settings_c, attr_name + "_interface", func_interface)


class SolverSettings(Settings):

    c_struct = SolverSettingsC
    init_c_struct = lib.init_solver_settings

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    logger: Optional[Callable[[str], None]]
    precond_interface: Any
    project_interface: Any
    conv_check_interface: Any
    logger_interface: Any


class StabilitySettings(Settings):

    c_struct = StabilitySettingsC
    init_c_struct = lib.init_stability_settings

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    logger: Optional[Callable[[str], None]]
    precond_interface: Any
    project_interface: Any
    logger_interface: Any


def auto_bind_fields(cls: type[Settings]):
    """
    this function automatically generates properties for non-pointer fields in a
    settings_c object that is an attribute of cls
    """
    for field_info in cls.c_struct._fields_:
        field_name, field_type = field_info[:2]

        # skip if function pointer, these will be initialized separately
        if field_type is c_void_p:
            continue

        # character arrays need to be handled separately
        elif field_type is c_char * (kw_len + 1):

            def make_property(name):

                def getter(self):
                    # get bytes from ctypes Structure
                    raw = getattr(self.settings_c, name)

                    # only get elements before null bytes and convert to strings
                    return raw.split(b"\0", 1)[0].decode("utf-8")

                def setter(self, value):
                    # check for type (will be a string if supplied by a user and a
                    # bytes when extracted from a ctypes Structure)
                    if isinstance(value, str):
                        b = value.encode("utf-8")[:kw_len]
                    elif isinstance(value, bytes):
                        b = value[:kw_len]

                    # pad with null bytes up to kw_len
                    b = b.ljust(kw_len, b"\0")

                    # copy to memory in ctypes Structure
                    offset = getattr(type(self.settings_c), name).offset
                    memmove(addressof(self.settings_c) + offset, b, kw_len)

                return property(getter, setter)

        else:

            def make_property(name):

                def getter(self):
                    return getattr(self.settings_c, name)

                def setter(self, value):
                    setattr(self.settings_c, name, value)

                return property(getter, setter)

        # capture field_name in current loop scope
        setattr(cls, field_name, make_property(field_name))


# ensure that appropriate fields are automatically set in settings_c object
auto_bind_fields(SolverSettings)
auto_bind_fields(StabilitySettings)


def solver(
    obj_func: Callable[[np.ndarray], float],
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    n_param: int,
    settings: SolverSettings,
):
    # define interfaces for callback functions
    @update_orbs_interface_type
    def update_orbs_interface(kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_funptr):
        """
        this function provides the interface to update the orbitals and to write the function value,
        gradient, Hessian diagonal to the memory provided through pointers and returns
        a function pointer to the Hessian linear transformation
        """
        # convert orbital rotation matrix pointer to numpy array
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(n_param,))
        grad = np.ctypeslib.as_array(grad_ptr, shape=(n_param,))
        h_diag = np.ctypeslib.as_array(h_diag_ptr, shape=(n_param,))

        # update orbitals and retrieve objective function, gradient, Hessian diagonal
        # and Hessian linear transformation function
        try:
            func_ptr[0], hess_x = update_orbs(kappa, grad, h_diag)
        except RuntimeError:
            return 1

        # attach the Hessian-vector product function to the solver function so that it
        # persists in Python to ensure that it is not garbage collected when the
        # current orbital updating function completes
        solver._hess_x_interface = hess_x_interface_factory(hess_x, n_param)

        # store the function pointer in hess_x_ptr so that it can be accessed by Fortran
        hess_x_funptr[0] = solver._hess_x_interface

        return 0

    @obj_func_interface_type
    def obj_func_interface(kappa_ptr, func_ptr):
        """
        this function interfaces the objective function
        """
        # convert orbital rotation matrix pointer to numpy array
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(n_param,))

        try:
            func_ptr[0] = obj_func(kappa)
        except RuntimeError:
            return 1

        return 0

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "precond", precond_interface_factory(settings.precond, n_param)
    )
    settings.set_optional_callback(
        "project", project_interface_factory(settings.project, n_param)
    )
    settings.set_optional_callback(
        "conv_check", conv_check_interface_factory(settings.conv_check)
    )
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    # define result and argument types
    lib.solver.restype = c_int
    lib.solver.argtypes = [
        update_orbs_interface_type,
        obj_func_interface_type,
        c_int,
        SolverSettingsC,
    ]

    # call Fortran function
    error = lib.solver(
        update_orbs_interface, obj_func_interface, n_param, settings.settings_c
    )

    if error:
        raise RuntimeError("OpenTrustRegion solver produced error.")


def stability_check(
    h_diag: np.ndarray,
    hess_x: Callable[[np.ndarray, np.ndarray], None],
    n_param: int,
    settings: StabilitySettings,
    kappa: Optional[np.ndarray] = None,
) -> bool:
    # define interfaces for callback functions
    hess_x_interface = hess_x_interface_factory(hess_x, n_param)

    # set interfaces for optional callback functions, these need to be set here since
    # the interface might need parameters that are not known when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "precond", precond_interface_factory(settings.precond, n_param)
    )
    settings.set_optional_callback(
        "project", project_interface_factory(settings.project, n_param)
    )
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    # define result and argument types
    lib.stability_check.restype = c_int
    lib.stability_check.argtypes = [
        POINTER(c_real),
        hess_x_interface_type,
        c_int,
        POINTER(c_bool),
        StabilitySettingsC,
        c_void_p,
    ]

    # initialize return variables
    stable = c_bool(False)

    # call Fortran function
    error = lib.stability_check(
        h_diag.ctypes.data_as(POINTER(c_real)),
        hess_x_interface,
        n_param,
        byref(stable),
        settings.settings_c,
        kappa.ctypes.data_as(POINTER(c_real)) if kappa is not None else kappa,
    )

    if error:
        raise RuntimeError("OpenTrustRegion stability check produced error.")

    return bool(stable)
