# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import unittest
from importlib import resources
from ctypes import CDLL, c_bool
from unittest.mock import patch

# check if numpy is available
try:
    import numpy as np

    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

# check if pyopentrustregion is installed or import module in same directory
try:
    from pyopentrustregion import solver, stability_check
except ImportError:
    from python_interface import solver, stability_check

# load the testsuite library
ext = "dylib" if sys.platform == "darwin" else "so"
try:
    with resources.path("pyopentrustregion", f"libtestsuite.{ext}") as lib_path:
        libtestsuite = CDLL(str(lib_path))
# fallback location if installation was not through setup.py
except OSError:
    try:
        fallback_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "../build", f"libtestsuite.{ext}")
        )
        libtestsuite = CDLL(fallback_path)
    except OSError:
        raise FileNotFoundError("Cannot find testsuite library.")


# define all tests
fortran_tests = {
    "opentrustregion_tests": [
        "truncated_conjugate_gradient",
        "level_shifted_davidson",
        "sanity_check",
        "accept_trust_region_step",
        "split_string_by_space",
        "log",
        "print_results",
        "minres",
        "jacobi_davidson_correction",
        "orthogonal_projection",
        "abs_diag_precond",
        "level_shifted_diag_precond",
        "set_default",
        "init_solver_settings",
        "init_stability_settings",
        "init_rng",
        "gram_schmidt",
        "generate_random_trial_vectors",
        "generate_trial_vectors",
        "add_column",
        "extend_symm_matrix",
        "min_eigval",
        "symm_mat_min_eig",
        "bracket",
        "newton_step",
        "bisection",
        "solver",
        "stability_check",
    ],
    "c_interface_tests": [
        "solver_c_wrapper",
        "stability_check_c_wrapper",
        "update_orbs_c_wrapper",
        "hess_x_c_wrapper",
        "obj_func_c_wrapper",
        "precond_c_wrapper",
        "conv_check_c_wrapper",
        "logger_c_wrapper",
        "set_default_c_ptr",
    ],
    "system_tests": ["h2o_atomic_fb", "h2o_saddle_fb"],
}

# define return type of Fortran functions
for tests in fortran_tests.values():
    for test in tests:
        getattr(libtestsuite, f"test_{test}").restype = c_bool


# define function to add tests to test classes
def add_tests(cls):

    def create_test(func_name):
        def test(self):
            result = getattr(libtestsuite, func_name)()
            if result:
                print(f" {func_name} PASSED")
            self.assertTrue(result, f"{func_name} failed")

        return test

    for func_name in cls.tests:
        setattr(cls, f"test_{func_name}", create_test(f"test_{func_name}"))

    return cls


@add_tests
class OpenTrustRegionTests(unittest.TestCase):
    """
    this class contains unit tests for opentrustregion
    """

    tests = fortran_tests["opentrustregion_tests"]

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running unit tests for OpenTrustRegion...")
        print(50 * "-")
        return super().setUpClass()


@add_tests
class CInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the C interface
    """

    tests = fortran_tests["c_interface_tests"]

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running unit tests for C interface...")
        print(50 * "-")
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class PyInterfaceTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running unit tests for Python interface...")
        print(50 * "-")
        return super().setUpClass()

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.libopentrustregion.solver",
        libtestsuite.mock_solver,
    )
    def test_solver_py_interface(self):
        """
        this function tests the solver python interface
        """

        def mock_obj_func(kappa):
            """
            this function is a mock function for the objective function
            """
            return np.sum(kappa)

        def mock_update_orbs(kappa):
            """
            this function is a mock function for the orbital update function
            """
            func = np.sum(kappa)
            grad = 2 * kappa
            h_diag = 3 * kappa

            def hess_x(x):
                return 4 * x

            return func, grad, h_diag, hess_x

        def mock_precond(residual, mu):
            """
            this function is a mock function for the preconditioner function
            """
            return mu * residual

        def mock_conv_check():
            """
            this function is a mock function for the convergence check function
            """
            return True

        def mock_logger(message):
            """
            this function is a mock function for the logging function
            """
            nonlocal test_logger
            if message == "test":
                test_logger = True
            return

        # call solver python interface without optional arguments
        solver(mock_obj_func, mock_update_orbs, 3)

        # initialize logging boolean
        test_logger = False

        # call solver python interface with optional arguments
        solver(
            mock_obj_func,
            mock_update_orbs,
            3,
            precond=mock_precond,
            conv_check=mock_conv_check,
            stability=False,
            line_search=True,
            davidson=False,
            jacobi_davidson=False,
            prefer_jacobi_davidson=True,
            conv_tol=1e-3,
            n_random_trial_vectors=5,
            start_trust_radius=0.2,
            n_macro=300,
            n_micro=200,
            global_red_factor=1e-2,
            local_red_factor=1e-3,
            seed=33,
            verbose=3,
            logger=mock_logger,
        )

        # check if logger was called correctly
        if not test_logger:
            print("test_solver_py_interface failed: Called logging function wrong.")

        self.assertTrue(
            libtestsuite.test_solver_result() or not test_logger,
            "test_solver_py_interface failed",
        )
        print("test_solver_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.libopentrustregion.stability_check",
        libtestsuite.mock_stability_check,
    )
    def test_stability_check_py_interface(self):
        """
        this function tests the stability check python interface
        """
        h_diag = np.full(3, 3.0, dtype=np.float64)

        def mock_hess_x(x):
            return 4 * x

        def mock_precond(residual, mu):
            """
            this function is a mock function for the preconditioner function
            """
            return mu * residual

        def mock_logger(message):
            """
            this function is a mock function for the logging function
            """
            nonlocal test_logger
            if message == "test":
                test_logger = True
            return

        # call stability check python interface without optional arguments
        stable, kappa = stability_check(h_diag, mock_hess_x, 3)

        msg = (
            "test_stability_check_py_interface failed: Returned stability boolean "
            "wrong."
        )
        if stable:
            print(msg)
        self.assertFalse(stable, msg)
        msg = "test_stability_check_py_interface failed: Returned direction wrong."
        if not np.allclose(kappa, np.full(3, 1.0, dtype=np.float64)):
            print(msg)
        self.assertTrue(np.allclose(kappa, np.full(3, 1.0, dtype=np.float64)), msg)

        # initialize logging boolean
        test_logger = False

        # call stability check python interface with optional arguments
        stable, kappa = stability_check(
            h_diag,
            mock_hess_x,
            3,
            precond=mock_precond,
            jacobi_davidson=False,
            conv_tol=1e-3,
            n_random_trial_vectors=3,
            n_iter=50,
            verbose=3,
            logger=mock_logger,
        )

        # check if logger was called correctly
        if not test_logger:
            print(
                "test_stability_check_py_interface failed: Called logging function "
                "wrong."
            )

        # check if returned variables are correct
        if stable:
            print(
                "test_stability_check_py_interface failed: Returned stability boolean "
                "wrong."
            )
        if not np.allclose(kappa, np.full(3, 1.0, dtype=np.float64)):
            print("test_stability_check_py_interface failed: Returned direction wrong.")

        self.assertTrue(
            libtestsuite.test_stability_check_result()
            or not test_logger
            or stable
            or not np.allclose(kappa, np.full(3, 1.0, dtype=np.float64)),
            "test_stability_check_py_interface failed",
        )
        print("test_stability_check_py_interface PASSED")


@add_tests
class SystemTests(unittest.TestCase):
    """
    this class contains system tests for opentrustregion
    """

    tests = fortran_tests["system_tests"]

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running system tests for OpenTrustRegion...")
        print(50 * "-")
        return super().setUpClass()


if __name__ == "__main__":
    unittest.main(verbosity=0)
