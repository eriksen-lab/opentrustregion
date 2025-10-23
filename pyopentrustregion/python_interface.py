# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import os
import sys
import numpy as np
from importlib import resources
from ctypes import (
    CDLL,
    CFUNCTYPE,
    POINTER,
    byref,
    string_at,
    cast,
    c_double,
    c_int64,
    c_int32,
    c_bool,
    c_void_p,
    c_char_p,
    Structure,
)
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any

# load the opentrustregion library
ext = "dylib" if sys.platform == "darwin" else "so"
try:
    with resources.path("pyopentrustregion", f"libopentrustregion.{ext}") as lib_path:
        libopentrustregion = CDLL(str(lib_path))
# fallback location if installation was not through setup.py
except OSError:
    try:
        fallback_path = os.path.abspath(
            os.path.join(
                os.path.dirname(__file__), "../build", f"libopentrustregion.{ext}"
            )
        )
        libopentrustregion = CDLL(fallback_path)
    except OSError:
        raise FileNotFoundError("Cannot find opentrustregion library.")

# determine integer size used in library
libopentrustregion.is_ilp64.restype = c_bool
if libopentrustregion.is_ilp64():
    c_ip = c_int64
else:
    c_ip = c_int32

# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant
# data
hess_x_interface_type = CFUNCTYPE(c_ip, POINTER(c_double), POINTER(c_double))
update_orbs_interface_type = CFUNCTYPE(
    c_ip,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(hess_x_interface_type),
)
obj_func_interface_type = CFUNCTYPE(c_ip, POINTER(c_double), POINTER(c_double))
precond_interface_type = CFUNCTYPE(
    c_ip, POINTER(c_double), POINTER(c_double), POINTER(c_double)
)
conv_check_interface_type = CFUNCTYPE(c_ip, POINTER(c_bool))
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
        ("conv_check", c_void_p),
        ("logger", c_void_p),
        ("stability", c_bool),
        ("line_search", c_bool),
        ("davidson", c_bool),
        ("jacobi_davidson", c_bool),
        ("prefer_jacobi_davidson", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_double),
        ("start_trust_radius", c_double),
        ("global_red_factor", c_double),
        ("local_red_factor", c_double),
        ("n_random_trial_vectors", c_ip),
        ("n_macro", c_ip),
        ("n_micro", c_ip),
        ("seed", c_ip),
        ("verbose", c_ip),
    ]


class StabilitySettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("logger", c_void_p),
        ("jacobi_davidson", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_double),
        ("n_random_trial_vectors", c_ip),
        ("n_iter", c_ip),
        ("verbose", c_ip),
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

        # initializes all fields from the C struct
        for field_info in self.settings_c._fields_:
            field_name = field_info[0]
            setattr(self, field_name, getattr(self.settings_c, field_name))

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
    init_c_struct = libopentrustregion.init_solver_settings

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    logger: Optional[Callable[[str], None]]
    precond_interface: Any
    conv_check_interface: Any
    logger_interface: Any


class StabilitySettings(Settings):

    c_struct = StabilitySettingsC
    init_c_struct = libopentrustregion.init_stability_settings

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    logger: Optional[Callable[[str], None]]
    precond_interface: Any
    logger_interface: Any


def auto_bind_fields(cls: type[Settings]):
    """
    this function automatically generates Python properties for non-pointer fields in a
    settings_c object that is an attribute of cls, ensuring synchronization between
    Python attributes and the C struct
    """
    for field_info in cls.c_struct._fields_:
        field_name, field_type = field_info[:2]

        # skip if function pointer
        if field_type == c_void_p:
            continue

        def make_property(name):
            # get attribute from _name
            def getter(self):
                return getattr(self, f"_{name}")

            # set attribute in _name and set attribute  in ctypes.Structure object
            def setter(self, value):
                setattr(self, f"_{name}", value)
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

        # attach the function to some object that persists in Fortran to ensure that it
        # is not garbage collected when the current function completes
        grad_ptr._hess_x_interface = hess_x_interface_factory(hess_x, n_param)

        # store the function pointer in hess_x_ptr so that it can be accessed by Fortran
        hess_x_funptr[0] = grad_ptr._hess_x_interface

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
    # the interface might need parameters that are not know when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "precond", precond_interface_factory(settings.precond, n_param)
    )
    settings.set_optional_callback(
        "conv_check", conv_check_interface_factory(settings.conv_check)
    )
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    # define result and argument types
    libopentrustregion.solver.restype = c_ip
    libopentrustregion.solver.argtypes = [
        update_orbs_interface_type,
        obj_func_interface_type,
        c_ip,
        SolverSettingsC,
    ]

    # call Fortran function
    error = libopentrustregion.solver(
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
    # the interface might need parameters that are not know when the attribute to
    # settings is set (e.g. n_param)
    settings.set_optional_callback(
        "precond", precond_interface_factory(settings.precond, n_param)
    )
    settings.set_optional_callback("logger", logger_interface_factory(settings.logger))

    # define result and argument types
    libopentrustregion.stability_check.restype = c_ip
    libopentrustregion.stability_check.argtypes = [
        POINTER(c_double),
        hess_x_interface_type,
        c_ip,
        POINTER(c_bool),
        StabilitySettingsC,
        c_void_p,
    ]

    # initialize return variables
    stable = c_bool(False)

    # call Fortran function
    error = libopentrustregion.stability_check(
        h_diag.ctypes.data_as(POINTER(c_double)),
        hess_x_interface,
        n_param,
        byref(stable),
        settings.settings_c,
        kappa.ctypes.data_as(POINTER(c_double)) if kappa is not None else kappa,
    )

    if error:
        raise RuntimeError("OpenTrustRegion stability check produced error.")

    return bool(stable)
