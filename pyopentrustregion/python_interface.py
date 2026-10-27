# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from __future__ import annotations

import os
import sys
import copy
import signal
import threading
import warnings
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
    c_size_t,
    sizeof,
    Structure,
)
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional, Any, Dict, Type, Union


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
conda_prefix = Path("/opt/anaconda1anaconda2anaconda3")


def _load_lib(lib_names: list) -> CDLL:
    """
    this function loads the first of the given libraries it finds, searching the
    conda prefix, the installed package, the directory of this file and ../build in
    this order
    """
    # try to load from installed package (conda)
    for lib_name in lib_names:
        conda_path_unix = conda_prefix / "lib" / lib_name
        conda_path_wind = conda_prefix / "Library" / "bin" / lib_name
        if conda_path_unix.exists():
            return CDLL(str(conda_path_unix))
        elif conda_path_wind.exists():
            return CDLL(str(conda_path_wind))

    # fallback: try to load from installed package (site-packages)
    for lib_name in lib_names:
        lib_path = resources.files("pyopentrustregion") / lib_name
        if lib_path.is_file():
            return CDLL(str(lib_path))

    # fallback: try to load from the directory of this file
    for lib_name in lib_names:
        local_path = os.path.join(os.path.dirname(__file__), lib_name)
        if os.path.exists(local_path):
            return CDLL(local_path)

    # fallback: try to load from ../build (development build)
    for lib_name in lib_names:
        build_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../build", lib_name)
        )
        if os.path.exists(build_path):
            return CDLL(build_path)

    # if all failed
    raise FileNotFoundError(f"Cannot find any of the libraries {lib_names}")


lib = _load_lib(lib_candidates)

# determine integer size used in library
ilp64 = c_bool.in_dll(lib, "ilp64")
c_int: Type[Union[c_int32, c_int64]]
if ilp64.value:
    c_int = c_int64
else:
    c_int = c_int32

# define real type
c_real = c_double

# fixed size strings for keywords
kw_len = c_int.in_dll(lib, "kw_len_c").value


def callback_type(restype, *argtypes):
    """
    this function builds the ctypes signature of a callback function, every callback
    function of the C interface receives the opaque host context as its last argument
    so it is appended here instead of at every declaration
    """
    return CFUNCTYPE(restype, *argtypes, c_void_p)


# callback function ctypes specifications, ctypes can only deal with simple return
# types so we interface to Fortran subroutines by creating pointers to the relevant data
hess_x_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
update_orbs_interface_type = callback_type(
    c_int,
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(c_real),
    POINTER(hess_x_interface_type),
)
obj_func_interface_type = callback_type(c_int, POINTER(c_real), POINTER(c_real))
precond_interface_type = callback_type(
    c_int, POINTER(c_real), POINTER(c_real), POINTER(c_real)
)
project_interface_type = callback_type(c_int, POINTER(c_real))
get_extra_trial_vectors_interface_type = callback_type(c_int, POINTER(c_real), c_int)
conv_check_interface_type = callback_type(c_int, POINTER(c_bool))
logger_interface_type = callback_type(None, c_char_p)


class InterruptionState(threading.local):
    """
    this class keeps, separately for every thread, a stack of the running calls, True
    for a callback function and False for the library, innermost last, and an
    interruption that arrived while the library ran and is raised later
    """

    stack: tuple = ()
    pending: Optional[BaseException] = None


interruption_state = InterruptionState()


def defer_interruption(signum, frame):
    """
    this function handles SIGINT during a call into the library, it raises a
    KeyboardInterrupt right away while a callback function runs, but only records it
    while the library runs, since Python would otherwise raise it on entering the next
    callback function, where ctypes would drop it
    """
    if interruption_state.stack and not interruption_state.stack[-1]:
        interruption_state.pending = KeyboardInterrupt()
    else:
        signal.default_int_handler(signum, frame)


# define interface factories
class CallbackInterface:
    """
    this class is the base class of the callback function interfaces, the C interface
    hands every callback function the opaque host context as its last argument while
    Python hosts carry their data in closures instead, so it is dropped here and the
    subclasses only implement call without it, an exception is reported to the library
    as an error, since ctypes would drop it
    """

    exception: Dict[str, BaseException]

    def __call__(self, *args):
        # catch every exception, since ctypes drops one that leaves this function and
        # hands the library an undefined return value instead
        outer = interruption_state.stack
        try:
            interruption_state.stack = outer + (True,)

            # stop the call for an interruption that arrived while the library ran
            if interruption_state.pending is not None:
                pending, interruption_state.pending = interruption_state.pending, None
                raise pending

            result = self.call(*args[:-1])
        except BaseException as e:
            # unmark the callback function before calling fail, so that an interruption
            # arriving during that call is recorded instead of raised outside the try
            interruption_state.stack = outer
            return self.fail(e)
        interruption_state.stack = outer
        return result

    def fail(self, e: BaseException) -> Optional[int]:
        """
        this function keeps the exception of a callback function, which is raised once
        the call has returned, and returns the error code that stops the call
        """
        self.exception["exc"] = e
        return 1


@dataclass
class ObjFuncInterface(CallbackInterface):
    """
    this class provides the interface for the objective function
    """

    obj_func: Callable[[np.ndarray], float]
    n_param: int
    exception: Dict[str, BaseException]

    def call(self, kappa_ptr, func_ptr) -> int:
        # convert matrix pointers to numpy arrays and make read-only
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(self.n_param,))
        kappa.flags.writeable = False

        func_ptr[0] = self.obj_func(kappa)

        return 0


@dataclass
class HessXInterface(CallbackInterface):
    """
    this class provides the interface for the Hessian linear transformation
    """

    hess_x: Callable[[np.ndarray, np.ndarray], None]
    n_param: int
    exception: Dict[str, BaseException]

    def call(self, x_ptr, hx_ptr) -> int:
        # convert trial vector pointer to numpy array and make read-only
        x = np.ctypeslib.as_array(x_ptr, shape=(self.n_param,))
        x.flags.writeable = False
        hx = np.ctypeslib.as_array(hx_ptr, shape=(self.n_param,))

        # perform linear transformation
        self.hess_x(x, hx)

        return 0


@dataclass
class UpdateOrbsInterface(CallbackInterface):
    """
    this class provides the interface to the orbital updating function
    """

    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ]
    n_param: int
    exception: Dict[str, BaseException]
    hess_x_funptr: Optional[Any] = None

    def call(self, kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_funptr) -> int:
        # convert matrix pointers to numpy arrays and make read-only
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(self.n_param,))
        kappa.flags.writeable = False
        grad = np.ctypeslib.as_array(grad_ptr, shape=(self.n_param,))
        h_diag = np.ctypeslib.as_array(h_diag_ptr, shape=(self.n_param,))

        # update orbitals and retrieve objective function, gradient, Hessian diagonal
        # and Hessian linear transformation function
        func_ptr[0], hess_x = self.update_orbs(kappa, grad, h_diag)

        # an orbital update which does not provide a Hessian linear transformation
        # leaves the pointer unset so that the solver reports the missing one, while
        # one that is not callable is reported as an error of the orbital update
        if hess_x is None:
            return 0
        if not callable(hess_x):
            self.exception["exc"] = TypeError(
                "Orbital update returned a Hessian linear transformation of type "
                f"{type(hess_x).__name__}, which is not callable."
            )
            return 1

        # keep the interface of the Hessian linear transformation on this object so
        # that it is not garbage collected while the library still calls it
        self.hess_x_funptr = hess_x_interface_type(
            HessXInterface(hess_x, self.n_param, self.exception)
        )
        hess_x_funptr[0] = self.hess_x_funptr

        return 0


@dataclass
class PrecondInterface(CallbackInterface):
    """
    this class provides the interface to the preconditioning function
    """

    precond: Callable[[np.ndarray, float, np.ndarray], None]
    n_param: int
    exception: Dict[str, BaseException]

    def call(self, residual_ptr, mu_ptr, precond_residual_ptr) -> int:
        # convert pointers to numpy arrays and float and make read-only
        residual = np.ctypeslib.as_array(residual_ptr, shape=(self.n_param,))
        residual.flags.writeable = False
        mu = mu_ptr[0]
        precond_residual = np.ctypeslib.as_array(
            precond_residual_ptr, shape=(self.n_param,)
        )

        # call preconditioner
        self.precond(residual, mu, precond_residual)

        return 0


@dataclass
class ProjectInterface(CallbackInterface):
    """
    this class provides the interface to the projection function
    """

    project: Callable[[np.ndarray], None]
    n_param: int
    exception: Dict[str, BaseException]

    def call(self, vector_ptr) -> int:
        # convert matrix pointers to numpy arrays
        vector = np.ctypeslib.as_array(vector_ptr, shape=(self.n_param,))

        # call projection function
        self.project(vector)

        return 0


@dataclass
class GetExtraTrialVectorsInterface(CallbackInterface):
    """
    this class provides the interface to the extra trial vector function
    """

    get_extra_trial_vectors: Callable[[np.ndarray], None]
    n_param: int
    exception: Dict[str, BaseException]

    def call(self, trial_vectors_ptr, n_extra_trial_vectors) -> int:
        # convert matrix pointer to numpy array with one trial vector per row
        trial_vectors = np.ctypeslib.as_array(
            trial_vectors_ptr, shape=(n_extra_trial_vectors, self.n_param)
        )

        # call extra trial vector function
        self.get_extra_trial_vectors(trial_vectors)

        return 0


@dataclass
class ConvCheckInterface(CallbackInterface):
    """
    this class provides the interface to the convergence check function
    """

    conv_check: Callable[[], bool]
    exception: Dict[str, BaseException]

    def call(self, conv_ptr) -> int:
        # call convergence check
        conv_ptr[0] = self.conv_check()

        return 0


@dataclass
class LoggerInterface(CallbackInterface):
    """
    this class provides the interface to the logging function
    """

    logger: Callable[[str], None]
    exception: Dict[str, BaseException]

    def call(self, message):
        # call logger
        self.logger(string_at(message).decode("utf-8"))

    def fail(self, e: BaseException) -> None:
        """
        this function keeps the first exception of the logging function, which returns
        no error code since the result does not depend on it, for a warning once the
        call has returned, and records an interruption, so that the next callback
        function stops the call
        """
        if isinstance(e, Exception):
            self.exception.setdefault("logger", e)
        else:
            interruption_state.pending = e


# define classes corresponding to C structs for settings
class StabilitySettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("project", c_void_p),
        ("get_extra_trial_vectors", c_void_p),
        ("conv_check", c_void_p),
        ("logger", c_void_p),
        ("hess_symm", c_bool),
        ("stop_on_instability", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_real),
        ("stability_thresh", c_real),
        ("n_random_trial_vectors", c_int),
        ("n_extra_trial_vectors", c_int),
        ("n_iter", c_int),
        ("jacobi_davidson_start", c_int),
        ("seed", c_int),
        ("verbose", c_int),
        ("diag_solver", c_char * (kw_len + 1)),
        ("context", c_void_p),
    ]


class SolverSettingsC(Structure):
    _fields_ = [
        ("precond", c_void_p),
        ("project", c_void_p),
        ("get_extra_trial_vectors", c_void_p),
        ("conv_check", c_void_p),
        ("logger", c_void_p),
        ("stability", c_bool),
        ("line_search", c_bool),
        ("grad_noise", c_bool),
        ("hess_symm", c_bool),
        ("initialized", c_bool),
        ("conv_tol", c_real),
        ("start_trust_radius", c_real),
        ("global_red_factor", c_real),
        ("local_red_factor", c_real),
        ("n_random_trial_vectors", c_int),
        ("n_extra_trial_vectors", c_int),
        ("n_macro", c_int),
        ("n_micro", c_int),
        ("jacobi_davidson_start", c_int),
        ("seed", c_int),
        ("verbose", c_int),
        ("subsystem_solver", c_char * (kw_len + 1)),
        ("stability_settings", StabilitySettingsC),
        ("context", c_void_p),
    ]


# define classes corresponding to C structs for states
class StabilityStateC(Structure):
    _fields_ = [
        ("min_eigval", c_real),
        ("residual_rms", c_real),
        ("iteration", c_int),
        ("n_hess_x", c_int),
    ]


class SolverStateC(Structure):
    _fields_ = [
        ("max_precision_reached", c_bool),
        ("func", c_real),
        ("grad_rms", c_real),
        ("trust_radius", c_real),
        ("relative_residual", c_real),
        ("level_shift", c_real),
        ("step_size", c_real),
        ("macro_iteration", c_int),
        ("micro_iteration", c_int),
        ("jacobi_davidson_iteration", c_int),
        ("n_update_orbs", c_int),
        ("n_obj_func", c_int),
        ("n_hess_x", c_int),
        ("stability_state", StabilityStateC),
    ]


def check_struct_sizes(library: CDLL) -> None:
    """
    this function checks that the settings and state structures have the sizes the
    library exports, so that a library of another version, such as a separately
    installed one or a stale build, is not read and written with the wrong layout,
    only a layout change that keeps the sizes goes undetected
    """
    for name, structure in [
        ("solver_settings", SolverSettingsC),
        ("stability_settings", StabilitySettingsC),
        ("solver_state", SolverStateC),
        ("stability_state", StabilityStateC),
    ]:
        # a library older than this check exports no sizes
        size_function = getattr(library, f"{name}_size", None)
        if size_function is not None:
            size_function.restype = c_size_t
        if size_function is None or size_function() != sizeof(structure):
            raise ImportError(
                f"The {name.replace('_', ' ')} of pyopentrustregion do not match those "
                f"of the loaded library {library._name}, which therefore comes from "
                "another version of OpenTrustRegion. Reinstall or rebuild the library "
                "together with pyopentrustregion."
            )


check_struct_sizes(lib)


# define classes for settings
class Settings:

    c_struct: type[Structure]
    init_c_struct: Any

    # attributes the settings classes set themselves besides the settings and the
    # optional callback functions
    own_attributes = ("settings_c", "_stability_settings")

    def __init__(self, settings_c: Optional[Structure] = None):
        """
        this function initializes the settings class
        """
        # specify C interface signature
        self.init_c_struct.argtypes = [POINTER(self.c_struct)]
        self.init_c_struct.restype = None

        # call C-side initialization to populate defaults, use passed settings if
        # provided
        if settings_c is None:
            self.settings_c = self.c_struct()
            self.init_c_struct(byref(self.settings_c))
        else:
            self.settings_c = settings_c

        # initializes all optional function pointers to None
        for name in self.callback_names():
            setattr(self, name, None)

    @classmethod
    def callback_names(cls) -> list:
        """
        this function returns the names of the optional callback functions, the
        function pointers of the C struct apart from the host context, which Python
        does not use
        """
        return [
            field_info[0]
            for field_info in cls.c_struct._fields_
            if field_info[1] is c_void_p and field_info[0] != "context"
        ]

    def __setattr__(self, name: str, value: Any) -> None:
        """
        this function sets an attribute only if it is a setting, an optional callback
        function or one the class sets itself, so that a misspelled setting raises
        instead of being silently ignored
        """
        if not (
            isinstance(getattr(type(self), name, None), property)
            or name in self.callback_names()
            or name in self.own_attributes
        ):
            raise AttributeError(f"{type(self).__name__} has no setting {name!r}.")
        super().__setattr__(name, value)

    def __copy__(self) -> Settings:
        """
        this function returns a copy of the settings with a C struct of its own, so
        that changing either leaves the other unchanged, which shares the callback
        functions
        """
        return self.copy_with(lambda func: func)

    def __deepcopy__(self, memo: dict) -> Settings:
        """
        this function returns a copy of the settings with a C struct of its own and
        deep copies of the callback functions
        """
        return self.copy_with(lambda func: copy.deepcopy(func, memo))

    def copy_with(self, copy_callback: Callable[[Any], Any]) -> Settings:
        """
        this function returns a copy of the settings with a C struct of its own, whose
        callback functions are copied by the given function
        """
        settings = type(self)()
        memmove(
            addressof(settings.settings_c),
            addressof(self.settings_c),
            sizeof(self.c_struct),
        )
        self.copy_callbacks(settings, copy_callback)
        return settings

    def copy_callbacks(
        self, settings: Settings, copy_callback: Callable[[Any], Any]
    ) -> None:
        """
        this function sets the optional callback functions of other settings to copies
        of these, made by the given function
        """
        for name in self.callback_names():
            setattr(settings, name, copy_callback(getattr(self, name)))

    def set_optional_callback(
        self,
        settings_c,
        interfaces,
        attr_name,
        func,
        func_interface,
        func_interface_type,
        *args,
    ):
        """
        this function sets a callback in the C struct of a call from a Python function
        and keeps the interface alive in the interfaces of the call
        """
        # create interface if function is not None
        interface = (
            func_interface_type(func_interface(func, *args))
            if func is not None
            else None
        )

        # create c_void_p pointer, points to NULL if func_interface is None
        setattr(settings_c, attr_name, cast(interface, c_void_p))

        # keep interface alive until the call returns
        interfaces.append(interface)


class SolverSettings(Settings):

    c_struct = SolverSettingsC
    init_c_struct = lib.init_solver_settings

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    get_extra_trial_vectors: Optional[Callable[[np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    logger: Optional[Callable[[str], None]]

    def __init__(self):
        super().__init__()
        self._stability_settings = StabilitySettings(
            settings_c=self.settings_c.stability_settings
        )

    @property
    def stability_settings(self) -> StabilitySettings:
        return self._stability_settings

    def copy_callbacks(
        self, settings: Settings, copy_callback: Callable[[Any], Any]
    ) -> None:
        """
        this function sets the optional callback functions of other settings, those of
        the nested stability check settings included, to copies of these, made by the
        given function
        """
        super().copy_callbacks(settings, copy_callback)
        self.stability_settings.copy_callbacks(
            settings.stability_settings, copy_callback
        )

    def set_optional_callbacks(
        self,
        settings_c: Structure,
        interfaces: list,
        n_param: int,
        exception: Dict[str, BaseException],
    ):
        """
        this function sets the interfaces for the optional callback functions in the C
        struct of a call
        """
        self.set_optional_callback(
            settings_c,
            interfaces,
            "precond",
            self.precond,
            PrecondInterface,
            precond_interface_type,
            n_param,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "project",
            self.project,
            ProjectInterface,
            project_interface_type,
            n_param,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "get_extra_trial_vectors",
            self.get_extra_trial_vectors,
            GetExtraTrialVectorsInterface,
            get_extra_trial_vectors_interface_type,
            n_param,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "conv_check",
            self.conv_check,
            ConvCheckInterface,
            conv_check_interface_type,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "logger",
            self.logger,
            LoggerInterface,
            logger_interface_type,
            exception,
        )
        self.stability_settings.set_optional_callbacks(
            settings_c.stability_settings, interfaces, n_param, exception
        )


class StabilitySettings(Settings):

    c_struct = StabilitySettingsC
    init_c_struct = lib.init_stability_settings

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    get_extra_trial_vectors: Optional[Callable[[np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    logger: Optional[Callable[[str], None]]

    def set_optional_callbacks(
        self,
        settings_c: Structure,
        interfaces: list,
        n_param: int,
        exception: Dict[str, BaseException],
    ):
        """
        this function sets the interfaces for the optional callback functions in the C
        struct of a call
        """
        self.set_optional_callback(
            settings_c,
            interfaces,
            "precond",
            self.precond,
            PrecondInterface,
            precond_interface_type,
            n_param,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "project",
            self.project,
            ProjectInterface,
            project_interface_type,
            n_param,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "get_extra_trial_vectors",
            self.get_extra_trial_vectors,
            GetExtraTrialVectorsInterface,
            get_extra_trial_vectors_interface_type,
            n_param,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "conv_check",
            self.conv_check,
            ConvCheckInterface,
            conv_check_interface_type,
            exception,
        )
        self.set_optional_callback(
            settings_c,
            interfaces,
            "logger",
            self.logger,
            LoggerInterface,
            logger_interface_type,
            exception,
        )


def auto_bind_fields(cls: type[Settings]):
    """
    this function automatically generates properties for non-pointer fields in a
    settings_c object that is an attribute of cls
    """
    for field_info in cls.c_struct._fields_:
        field_name, field_type = field_info[:2]

        # skip if function pointer (these will be initialized separately), nested
        # structures or the flag of initialized C settings, which Python settings
        # always are
        if (
            field_type is c_void_p
            or (isinstance(field_type, type) and issubclass(field_type, Structure))
            or field_name == "initialized"
        ):
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


# define classes for states, which the library keeps current during a call and which
# can therefore be read by the callback functions
class State:

    c_struct: type[Structure]

    def __init__(self, state_c: Optional[Structure] = None):
        """
        this function initializes the state class on a new C struct, which holds the
        default values, or on a given one
        """
        self.state_c = self.c_struct() if state_c is None else state_c

    def __copy__(self) -> State:
        """
        this function returns a copy of the state with a C struct of its own, so that
        a later call writing into the state leaves the copy unchanged
        """
        state = type(self)()
        memmove(
            addressof(state.state_c), addressof(self.state_c), sizeof(self.c_struct)
        )
        return state

    def __deepcopy__(self, memo: dict) -> State:
        """
        this function returns a copy of the state with a C struct of its own, a state
        holds nothing else that could be copied
        """
        return self.__copy__()


class SolverState(State):

    c_struct = SolverStateC

    def __init__(self):
        super().__init__()
        self._stability_state = StabilityState(state_c=self.state_c.stability_state)

    @property
    def stability_state(self) -> StabilityState:
        return self._stability_state


class StabilityState(State):

    c_struct = StabilityStateC


def bind_state_fields(cls: type[State]):
    """
    this function generates read-only properties for the fields of the C struct of a
    state class, apart from the nested state
    """
    for field_info in cls.c_struct._fields_:
        field_name, field_type = field_info[:2]

        # skip nested structures
        if isinstance(field_type, type) and issubclass(field_type, Structure):
            continue

        def make_property(name):

            def getter(self):
                return getattr(self.state_c, name)

            return property(getter)

        # capture field_name in current loop scope
        setattr(cls, field_name, make_property(field_name))


# ensure that the fields of the state_c object can be read
bind_state_fields(SolverState)
bind_state_fields(StabilityState)


def raise_on_failure(error: int, exception: Dict[str, BaseException], name: str):
    """
    this function raises if a call into the library failed, with the exception raised
    by a failing callback function as its cause, a failed logger cannot affect the
    result and is therefore only reported as a warning, both before an error and
    instead of one, while an interruption of a callback function, such as a
    KeyboardInterrupt, is raised as it is
    """
    if "exc" in exception and not isinstance(exception["exc"], Exception):
        raise exception["exc"]
    if "logger" in exception:
        warning = RuntimeWarning(
            f"OpenTrustRegion {name} called a logging function that raised "
            f"{exception['logger']!r}, so log messages from that point on may be "
            "missing."
        )
        warning.__cause__ = exception["logger"]
        warnings.warn(warning, stacklevel=3)
    if error:
        raise RuntimeError(
            f"OpenTrustRegion {name} produced error (code {error})."
        ) from exception.get("exc")


def call_library(function: Callable[..., int], *args) -> int:
    """
    this function calls into the library and marks it as innermost while it runs, in
    the main thread it also replaces Python's default handler of SIGINT with
    defer_interruption meanwhile, and raises an interruption that arrived after the
    last callback function once the call has returned
    """
    outer = interruption_state.stack
    install = (
        threading.current_thread() is threading.main_thread()
        and signal.getsignal(signal.SIGINT) is signal.default_int_handler
    )
    try:
        if install:
            signal.signal(signal.SIGINT, defer_interruption)
        interruption_state.stack = outer + (False,)
        error = function(*args)
    finally:
        # restore the handler while the library is still marked as innermost, since
        # Python handles a pending signal while it switches the handler, which
        # defer_interruption then records instead of raising it here
        if install:
            signal.signal(signal.SIGINT, signal.default_int_handler)
        interruption_state.stack = outer
    if interruption_state.pending is not None:
        pending, interruption_state.pending = interruption_state.pending, None
        raise pending
    return error


def solver(
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    obj_func: Callable[[np.ndarray], float],
    n_param: int,
    settings: Optional[SolverSettings] = None,
    state: Optional[SolverState] = None,
):
    # variable to capture exceptions in the callback interfaces
    exception: Dict[str, BaseException] = {}

    # define interfaces for callback functions
    obj_func_interface = obj_func_interface_type(
        ObjFuncInterface(obj_func, n_param, exception)
    )
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, n_param, exception)
    )

    # set interfaces for optional callback functions in a copy of the C struct of the
    # settings, as the C entry points do, so that a call never writes into the settings
    # object and calls running at the same time can share it, and keep them alive until
    # the call returns, these need to be set here since the interface might need
    # parameters that are not known when the attribute to settings is set (e.g. n_param)
    settings_c = None
    interfaces: list = []
    if settings is not None:
        settings_c = SolverSettingsC.from_buffer_copy(settings.settings_c)
        settings.set_optional_callbacks(settings_c, interfaces, n_param, exception)

    # define result and argument types
    lib.solver.restype = c_int
    lib.solver.argtypes = [
        update_orbs_interface_type,
        obj_func_interface_type,
        c_int,
        POINTER(SolverSettingsC),
        POINTER(SolverStateC),
    ]

    # call Fortran function
    error = call_library(
        lib.solver,
        update_orbs_interface,
        obj_func_interface,
        n_param,
        byref(settings_c) if settings_c is not None else None,
        byref(state.state_c) if state is not None else None,
    )

    raise_on_failure(error, exception, "solver")


def stability_check(
    h_diag: np.ndarray,
    hess_x: Callable[[np.ndarray, np.ndarray], None],
    n_param: int,
    settings: Optional[StabilitySettings] = None,
    state: Optional[StabilityState] = None,
    kappa: Optional[np.ndarray] = None,
) -> bool:
    # convert the Hessian diagonal into a contiguous array of reals, check its shape
    # and ensure the lowest eigenvector array is writable, contiguous and has the right
    # type and shape
    h_diag = np.ascontiguousarray(h_diag, dtype=np.float64)
    if h_diag.shape != (n_param,):
        raise ValueError(
            f"Hessian diagonal has shape {h_diag.shape} instead of ({n_param},)."
        )
    if kappa is not None and not (
        isinstance(kappa, np.ndarray)
        and kappa.dtype == np.float64
        and kappa.shape == (n_param,)
        and kappa.flags.c_contiguous
        and kappa.flags.writeable
    ):
        raise ValueError(
            "The array for the lowest eigenvector has to be a writeable, contiguous "
            f"float64 array of shape ({n_param},)."
        )

    # variable to capture exceptions in the callback interfaces
    exception: Dict[str, BaseException] = {}

    # define interfaces for callback functions
    hess_x_interface = hess_x_interface_type(HessXInterface(hess_x, n_param, exception))

    # set interfaces for optional callback functions in a copy of the C struct of the
    # settings, as the C entry points do, so that a call never writes into the settings
    # object and calls running at the same time can share it, and keep them alive until
    # the call returns, these need to be set here since the interface might need
    # parameters that are not known when the attribute to settings is set (e.g. n_param)
    settings_c = None
    interfaces: list = []
    if settings is not None:
        settings_c = StabilitySettingsC.from_buffer_copy(settings.settings_c)
        settings.set_optional_callbacks(settings_c, interfaces, n_param, exception)

    # define result and argument types
    lib.stability_check.restype = c_int
    lib.stability_check.argtypes = [
        POINTER(c_real),
        hess_x_interface_type,
        c_int,
        POINTER(c_bool),
        POINTER(StabilitySettingsC),
        POINTER(StabilityStateC),
        POINTER(c_real),
    ]

    # initialize return variable
    stable = c_bool(False)

    # call Fortran function
    error = call_library(
        lib.stability_check,
        h_diag.ctypes.data_as(POINTER(c_real)),
        hess_x_interface,
        n_param,
        byref(stable),
        byref(settings_c) if settings_c is not None else None,
        byref(state.state_c) if state is not None else None,
        kappa.ctypes.data_as(POINTER(c_real)) if kappa is not None else None,
    )

    raise_on_failure(error, exception, "stability check")

    return stable.value
