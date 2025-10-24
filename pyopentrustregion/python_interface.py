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
    c_double,
    c_long,
    c_bool,
    c_void_p,
    c_char_p,
)
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Tuple, Callable, Optional


# load the opentrustregion library, fallback to testsuite in case opentrustregion was
# statically compiled
ext = "dylib" if sys.platform == "darwin" else "so"
lib_candidates = [
    f"libopentrustregion.{ext}",
    f"libopentrustregion32.{ext}",
    f"libopentrustregion64.{ext}",
    f"libtestsuite.{ext}",
]
for lib_name in lib_candidates:
    try:
        with resources.path("pyopentrustregion", lib_name) as lib_path:
            lib = CDLL(str(lib_path))
            break
    except OSError:
        # fallback for non-installed or dev build
        try:
            lib = CDLL(
                os.path.abspath(
                    os.path.join(os.path.dirname(__file__), "../build", lib_name)
                )
            )
            break
        except OSError:
            pass
else:
    raise FileNotFoundError(
        f"Cannot find either opentrustregion or testsuite library ({lib_candidates})"
    )


def solver(
    obj_func: Callable[[np.ndarray], float],
    update_orbs: Callable[
        [np.ndarray, np.ndarray, np.ndarray],
        Tuple[float, Callable[[np.ndarray, np.ndarray], None]],
    ],
    n_param: int,
    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]] = None,
    conv_check: Optional[Callable[[None], bool]] = None,
    stability: Optional[bool] = None,
    line_search: Optional[bool] = None,
    davidson: Optional[bool] = None,
    jacobi_davidson: Optional[bool] = None,
    prefer_jacobi_davidson: Optional[bool] = None,
    conv_tol: Optional[float] = None,
    n_random_trial_vectors: Optional[int] = None,
    start_trust_radius: Optional[float] = None,
    n_macro: Optional[int] = None,
    n_micro: Optional[int] = None,
    global_red_factor: Optional[float] = None,
    local_red_factor: Optional[float] = None,
    seed: Optional[int] = None,
    verbose: Optional[int] = None,
    logger: Optional[Callable[[str], None]] = None,
):
    # callback function ctypes specifications, ctypes can only deal with simple return
    # types so we interface to Fortran subroutines by creating pointers to the relevant
    # data
    hess_x_interface_type = CFUNCTYPE(c_long, POINTER(c_double), POINTER(c_double))
    update_orbs_interface_type = CFUNCTYPE(
        c_long,
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(hess_x_interface_type),
    )
    obj_func_interface_type = CFUNCTYPE(c_long, POINTER(c_double), POINTER(c_double))
    precond_interface_type = CFUNCTYPE(
        c_long, POINTER(c_double), POINTER(c_double), POINTER(c_double)
    )
    conv_check_interface_type = CFUNCTYPE(c_long, POINTER(c_bool))
    logger_interface_type = CFUNCTYPE(None, c_char_p)

    @update_orbs_interface_type
    def update_orbs_interface(kappa_ptr, func_ptr, grad_ptr, h_diag_ptr, hess_x_funptr):
        """
        this function updates the orbitals and writes pointers to the function value,
        gradient, Hessian diagonal and returns a procedure function to the Hessian
        linear transformation
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

        @hess_x_interface_type
        def hess_x_interface(x_ptr, hx_ptr):
            """
            this function is returned by the orbital update function and calculates the
            Hessian linear transformation for a trial vector for the updated orbitals
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

        # attach the function to some object that persists in Fortran to ensure that it
        # is not garbage collected when the current function completes
        grad_ptr._hess_x_interface = hess_x_interface

        # store the function pointer in hess_x_ptr so that it can be accessed by Fortran
        hess_x_funptr[0] = hess_x_interface

        return 0

    @obj_func_interface_type
    def obj_func_interface(kappa_ptr, func_ptr):
        """
        this function returns the function value
        """
        # convert orbital rotation matrix pointer to numpy array
        kappa = np.ctypeslib.as_array(kappa_ptr, shape=(n_param,))

        try:
            func_ptr[0] = obj_func(kappa)
        except RuntimeError:
            return 1

        return 0

    @precond_interface_type
    def precond_interface(residual_ptr, mu_ptr, precond_residual_ptr):
        """
        this function returns the preconditioner
        """
        # convert pointers to numpy array and float
        residual = np.ctypeslib.as_array(residual_ptr, shape=(n_param,))
        mu = mu_ptr[0]
        precond_residual = np.ctypeslib.as_array(precond_residual_ptr, shape=(n_param,))

        # get preconditioner
        try:
            precond(residual, mu, precond_residual)
        except RuntimeError:
            return 1

        return 0

    @conv_check_interface_type
    def conv_check_interface(conv_ptr):
        """
        this function performs a convergence check
        """
        try:
            conv_ptr[0] = conv_check()
        except RuntimeError:
            return 1

        return 0

    @logger_interface_type
    def logger_interface(message):
        """
        this function logs results
        """
        # call logger
        logger(string_at(message).decode("utf-8"))

    # define result and argument types
    lib.solver.restype = c_long
    lib.solver.argtypes = [
        update_orbs_interface_type,
        obj_func_interface_type,
        c_long,
        c_void_p if precond is None else precond_interface_type,
        c_void_p if conv_check is None else conv_check_interface_type,
        POINTER(c_bool),
        POINTER(c_bool),
        POINTER(c_bool),
        POINTER(c_bool),
        POINTER(c_bool),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_long),
        POINTER(c_double),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_long),
        c_void_p if logger is None else logger_interface_type,
    ]

    # call Fortran function
    error = lib.solver(
        update_orbs_interface,
        obj_func_interface,
        n_param,
        None if precond is None else precond_interface,
        None if conv_check is None else conv_check_interface,
        None if stability is None else byref(c_bool(stability)),
        None if line_search is None else byref(c_bool(line_search)),
        None if davidson is None else byref(c_bool(davidson)),
        None if jacobi_davidson is None else byref(c_bool(jacobi_davidson)),
        (
            None
            if prefer_jacobi_davidson is None
            else byref(c_bool(prefer_jacobi_davidson))
        ),
        None if conv_tol is None else byref(c_double(conv_tol)),
        (
            None
            if n_random_trial_vectors is None
            else byref(c_long(n_random_trial_vectors))
        ),
        None if start_trust_radius is None else byref(c_double(start_trust_radius)),
        None if n_macro is None else byref(c_long(n_macro)),
        None if n_micro is None else byref(c_long(n_micro)),
        None if global_red_factor is None else byref(c_double(global_red_factor)),
        None if local_red_factor is None else byref(c_double(local_red_factor)),
        None if seed is None else byref(c_long(seed)),
        None if verbose is None else byref(c_long(verbose)),
        None if logger is None else logger_interface,
    )

    if error:
        raise RuntimeError("OpenTrustRegion solver produced error.")


def stability_check(
    h_diag: np.ndarray,
    hess_x: Callable[[np.ndarray, np.ndarray], None],
    n_param: int,
    precond: Optional[Callable[[np.ndarray, float, np.ndarray], None]] = None,
    jacobi_davidson: Optional[bool] = None,
    conv_tol: Optional[float] = None,
    n_random_trial_vectors: Optional[int] = None,
    n_iter: Optional[int] = None,
    verbose: Optional[int] = None,
    logger: Optional[Callable[[str], None]] = None,
) -> Tuple[bool, np.ndarray]:
    # callback function ctypes specifications, ctypes can only deal with simple return
    # types so we interface to Fortran subroutines by creating data to the relevant data
    hess_x_interface_type = CFUNCTYPE(c_long, POINTER(c_double), POINTER(c_double))
    precond_interface_type = CFUNCTYPE(
        c_long, POINTER(c_double), POINTER(c_double), POINTER(c_double)
    )
    logger_interface_type = CFUNCTYPE(None, c_char_p)

    @hess_x_interface_type
    def hess_x_interface(x_ptr, hx_ptr):
        # convert trial vector pointer to numpy array
        x = np.ctypeslib.as_array(x_ptr, shape=(n_param,))
        hx = np.ctypeslib.as_array(hx_ptr, shape=(n_param,))

        # perform linear transformation
        try:
            hess_x(x, hx)
        except RuntimeError:
            return 1

        return 0

    @precond_interface_type
    def precond_interface(residual_ptr, mu_ptr, precond_residual_ptr):
        # convert pointers to numpy array and float
        residual = np.ctypeslib.as_array(residual_ptr, shape=(n_param,))
        mu = mu_ptr[0]
        precond_residual = np.ctypeslib.as_array(precond_residual_ptr, shape=(n_param,))

        # perform linear transformation
        try:
            precond(residual, mu, precond_residual)
        except RuntimeError:
            return 1

        return 0

    @logger_interface_type
    def logger_interface(message):
        # call logger
        logger(string_at(message).decode("utf-8"))

    # define result and argument types
    lib.stability_check.restype = c_long
    lib.stability_check.argtypes = [
        POINTER(c_double),
        hess_x_interface_type,
        c_long,
        POINTER(c_bool),
        POINTER(c_double),
        c_void_p if precond is None else precond_interface_type,
        POINTER(c_bool),
        POINTER(c_double),
        POINTER(c_long),
        POINTER(c_long),
        POINTER(c_long),
        c_void_p if logger is None else logger_interface_type,
    ]

    # initialize return variables
    stable = c_bool(False)
    kappa = np.empty(n_param, dtype=np.float64)

    # call Fortran function
    error = lib.stability_check(
        h_diag.ctypes.data_as(POINTER(c_double)),
        hess_x_interface,
        n_param,
        byref(stable),
        kappa.ctypes.data_as(POINTER(c_double)),
        None if precond is None else precond_interface,
        None if jacobi_davidson is None else byref(c_bool(jacobi_davidson)),
        None if conv_tol is None else byref(c_double(conv_tol)),
        (
            None
            if n_random_trial_vectors is None
            else byref(c_long(n_random_trial_vectors))
        ),
        None if n_iter is None else byref(c_long(n_iter)),
        None if verbose is None else byref(c_long(verbose)),
        None if logger is None else logger_interface,
    )

    if error:
        raise RuntimeError("OpenTrustRegion stability check produced error.")

    return bool(stable), kappa
