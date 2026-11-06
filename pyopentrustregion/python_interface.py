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
from contextvars import ContextVar
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
    f"libotrtestsuite.{ext}",
]
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


# the call whose library runs innermost in the current thread, None while a callback
# function or no call runs
running_call: ContextVar[Optional[LibraryCall]] = ContextVar(
    "running_call", default=None
)


def defer_interruption(signum, frame):
    """
    this function handles SIGINT during a call into the library, it raises a
    KeyboardInterrupt right away while a callback function runs, but only records it
    with the call while the library runs, since Python would otherwise raise it on
    entering the next callback function, where ctypes would drop it
    """
    library_call = running_call.get()
    if library_call is not None:
        library_call.interruption = KeyboardInterrupt()
    else:
        signal.default_int_handler(signum, frame)


@dataclass(eq=False)
class LibraryCall:
    """
    this class makes a call into the library with a number of parameters, the
    interfaces of its callback functions keep what these raise here, since ctypes would
    drop it, the exception of a failed callback function, the first one of the logging
    function and an interruption, such as a KeyboardInterrupt of a callback function or
    SIGINT while the library runs, which the call raises once the library has returned
    """

    name: str
    n_param: int
    failure: Optional[Exception] = None
    logger_failure: Optional[Exception] = None
    interruption: Optional[BaseException] = None

    def run(self, function: Callable[..., int], *args) -> None:
        """
        this function calls into the library and marks it as innermost while it runs,
        in the main thread it also replaces Python's default handler of SIGINT with
        defer_interruption meanwhile, and reports what the callback functions raised
        once the library has returned
        """
        outer = running_call.get()
        install = (
            threading.current_thread() is threading.main_thread()
            and signal.getsignal(signal.SIGINT) is signal.default_int_handler
        )
        try:
            try:
                if install:
                    signal.signal(signal.SIGINT, defer_interruption)
                running_call.set(self)
                error = function(*args)
            finally:
                # restore the handler while the library is still marked as innermost,
                # since Python handles a pending signal while it switches the handler,
                # which defer_interruption then records instead of raising it here
                if install:
                    signal.signal(signal.SIGINT, signal.default_int_handler)
                running_call.set(outer)
            self.report(error)
        finally:
            # forget what the callback functions raised, since the interfaces in its
            # tracebacks refer to this call, which would otherwise form a reference
            # cycle that keeps the frames of the call and the callback functions alive
            # until the next garbage collection
            self.failure = self.logger_failure = self.interruption = None

    def report(self, error: int) -> None:
        """
        this function raises an interruption as it is and an error with the exception
        of a failed callback function as its cause, while a failed logging function
        cannot affect the result and is therefore only reported as a warning, both
        before an error and instead of one, in a frame of its own, since the tracebacks
        of the callback functions refer to the frame that called the library, which
        would otherwise form a reference cycle with the warning
        """
        if self.interruption is not None:
            raise self.interruption
        if self.logger_failure is not None:
            warning = RuntimeWarning(
                f"OpenTrustRegion {self.name} called a logging function that raised "
                f"{self.logger_failure!r}, so log messages from that point on may be "
                "missing."
            )
            warning.__cause__ = self.logger_failure
            warnings.warn(warning, stacklevel=4)
        if error:
            raise RuntimeError(
                f"OpenTrustRegion {self.name} produced error (code {error})."
            ) from self.failure


# define interface factories
class CallbackInterface:
    """
    this class is the base class of the callback function interfaces, the C interface
    hands every callback function the opaque host context as its last argument while
    Python hosts carry their data in closures instead, so it is dropped here and the
    subclasses only implement call without it, an exception is reported to the library
    as an error, since ctypes would drop it
    """

    library_call: LibraryCall

    @property
    def n_param(self) -> int:
        return self.library_call.n_param

    def __call__(self, *args):
        # catch every exception, since ctypes drops one that leaves this function and
        # hands the library an undefined return value instead, and mark the callback
        # function as innermost while it runs, so that SIGINT raises right away
        outer = running_call.get()
        try:
            running_call.set(None)

            # report an interruption that arrived while the library ran without calling
            # the callback function, the call raises it once the library has returned
            if self.library_call.interruption is not None:
                raise self.library_call.interruption

            result = self.call(*args[:-1])
        except BaseException as e:
            # unmark the callback function before calling fail, so that an interruption
            # arriving during that call is recorded instead of raised outside the try
            running_call.set(outer)
            return self.fail(e)
        running_call.set(outer)
        return result

    def fail(self, e: BaseException) -> Optional[int]:
        """
        this function keeps the exception of a callback function with the call, which
        raises it once the library has returned, and returns the error code that stops
        the library
        """
        if isinstance(e, Exception):
            self.library_call.failure = e
        else:
            self.library_call.interruption = e
        return 1


@dataclass
class ObjFuncInterface(CallbackInterface):
    """
    this class provides the interface for the objective function
    """

    obj_func: Callable[[np.ndarray], float]
    library_call: LibraryCall

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
    library_call: LibraryCall

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
    library_call: LibraryCall
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
            self.library_call.failure = TypeError(
                "Orbital update returned a Hessian linear transformation of type "
                f"{type(hess_x).__name__}, which is not callable."
            )
            return 1

        # keep the interface of the Hessian linear transformation on this object so
        # that it is not garbage collected while the library still calls it
        self.hess_x_funptr = hess_x_interface_type(
            HessXInterface(hess_x, self.library_call)
        )
        hess_x_funptr[0] = self.hess_x_funptr

        return 0


@dataclass
class PrecondInterface(CallbackInterface):
    """
    this class provides the interface to the preconditioning function
    """

    precond: Callable[[np.ndarray, float, np.ndarray], None]
    library_call: LibraryCall

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
    library_call: LibraryCall

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
    library_call: LibraryCall

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
    library_call: LibraryCall

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
    library_call: LibraryCall

    def call(self, message):
        # call logger
        self.logger(string_at(message).decode("utf-8"))

    def fail(self, e: BaseException) -> None:
        """
        this function keeps the first exception of the logging function, which returns
        no error code since the result does not depend on it, for a warning once the
        library has returned, and an interruption, so that the next callback function
        stops the library
        """
        if not isinstance(e, Exception):
            self.library_call.interruption = e
        elif self.library_call.logger_failure is None:
            self.library_call.logger_failure = e


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

    # the interface and the C function type of every callback function of the settings
    callback_fields: Dict[str, Tuple[type, Any]] = {}

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
        for name in self.callback_fields:
            setattr(self, name, None)

    def __setattr__(self, name: str, value: Any) -> None:
        """
        this function sets an attribute only if it is a setting, an optional callback
        function or one the class sets itself, so that a misspelled setting raises
        instead of being silently ignored
        """
        if not (
            isinstance(getattr(type(self), name, None), property)
            or name in self.callback_fields
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
        for name in self.callback_fields:
            setattr(settings, name, copy_callback(getattr(self, name)))

    def set_optional_callbacks(
        self, settings_c: Structure, interfaces: list, library_call: LibraryCall
    ):
        """
        this function sets the interfaces of the optional callback functions in the C
        struct of a call and keeps them alive in the interfaces of the call
        """
        for name, (interface, interface_type) in self.callback_fields.items():
            func = getattr(self, name)
            callback = (
                interface_type(interface(func, library_call))
                if func is not None
                else None
            )

            # write the address of the interface, NULL without one; ctypes.cast would
            # put the interface in a reference cycle, which keeps it and the function
            # it wraps alive after the call until the next garbage collection
            setattr(
                settings_c,
                name,
                c_void_p.from_buffer(callback).value if callback is not None else None,
            )
            interfaces.append(callback)


class SolverSettings(Settings):

    c_struct = SolverSettingsC
    init_c_struct = lib.init_solver_settings
    callback_fields = {
        "precond": (PrecondInterface, precond_interface_type),
        "project": (ProjectInterface, project_interface_type),
        "get_extra_trial_vectors": (
            GetExtraTrialVectorsInterface,
            get_extra_trial_vectors_interface_type,
        ),
        "conv_check": (ConvCheckInterface, conv_check_interface_type),
        "logger": (LoggerInterface, logger_interface_type),
    }

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
        self, settings_c: Structure, interfaces: list, library_call: LibraryCall
    ):
        """
        this function sets the interfaces of the optional callback functions in the C
        struct of a call, those of the nested stability check settings included
        """
        super().set_optional_callbacks(settings_c, interfaces, library_call)
        self.stability_settings.set_optional_callbacks(
            settings_c.stability_settings, interfaces, library_call
        )


class StabilitySettings(Settings):

    c_struct = StabilitySettingsC
    init_c_struct = lib.init_stability_settings
    callback_fields = {
        "precond": (PrecondInterface, precond_interface_type),
        "project": (ProjectInterface, project_interface_type),
        "get_extra_trial_vectors": (
            GetExtraTrialVectorsInterface,
            get_extra_trial_vectors_interface_type,
        ),
        "conv_check": (ConvCheckInterface, conv_check_interface_type),
        "logger": (LoggerInterface, logger_interface_type),
    }

    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]]
    project: Optional[Callable[[np.ndarray], None]]
    get_extra_trial_vectors: Optional[Callable[[np.ndarray], None]]
    conv_check: Optional[Callable[[], bool]]
    logger: Optional[Callable[[str], None]]


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
    # call into the library, which keeps what the callback functions raise
    library_call = LibraryCall("solver", n_param)

    # define interfaces for callback functions
    obj_func_interface = obj_func_interface_type(
        ObjFuncInterface(obj_func, library_call)
    )
    update_orbs_interface = update_orbs_interface_type(
        UpdateOrbsInterface(update_orbs, library_call)
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
        settings.set_optional_callbacks(settings_c, interfaces, library_call)

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
    library_call.run(
        lib.solver,
        update_orbs_interface,
        obj_func_interface,
        n_param,
        byref(settings_c) if settings_c is not None else None,
        byref(state.state_c) if state is not None else None,
    )


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

    # call into the library, which keeps what the callback functions raise
    library_call = LibraryCall("stability check", n_param)

    # define interfaces for callback functions
    hess_x_interface = hess_x_interface_type(HessXInterface(hess_x, library_call))

    # set interfaces for optional callback functions in a copy of the C struct of the
    # settings, as the C entry points do, so that a call never writes into the settings
    # object and calls running at the same time can share it, and keep them alive until
    # the call returns, these need to be set here since the interface might need
    # parameters that are not known when the attribute to settings is set (e.g. n_param)
    settings_c = None
    interfaces: list = []
    if settings is not None:
        settings_c = StabilitySettingsC.from_buffer_copy(settings.settings_c)
        settings.set_optional_callbacks(settings_c, interfaces, library_call)

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
    library_call.run(
        lib.stability_check,
        h_diag.ctypes.data_as(POINTER(c_real)),
        hess_x_interface,
        n_param,
        byref(stable),
        byref(settings_c) if settings_c is not None else None,
        byref(state.state_c) if state is not None else None,
        kappa.ctypes.data_as(POINTER(c_real)) if kappa is not None else None,
    )

    return stable.value
