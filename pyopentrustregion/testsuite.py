# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import unittest
from importlib import resources
from ctypes import (
    CDLL,
    c_bool,
    c_char,
    c_void_p,
    Array,
    byref,
    CFUNCTYPE,
    Structure,
    POINTER,
)
from unittest.mock import patch
from pathlib import Path

# check if numpy is available
try:
    import numpy as np

    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

# check if pyopentrustregion is installed or import module in same directory
try:
    from pyopentrustregion.python_interface import (
        SolverSettings,
        StabilitySettings,
        solver,
        stability_check,
        c_int,
        c_real,
        ext,
        conda_prefix,
    )
except ImportError:
    sys.path.insert(0, str(Path(__file__).parent.absolute()))
    from pyopentrustregion.python_interface import (
        SolverSettings,
        StabilitySettings,
        solver,
        stability_check,
        c_int,
        c_real,
        ext,
        conda_prefix,
    )

lib = None

# try to load from installed package (conda)
lib_name = f"libotrtestsuite.{ext}"
conda_path_unix = conda_prefix / "lib" / lib_name
conda_path_wind = conda_prefix / "Library" / "bin" / lib_name
if conda_path_unix.exists():
    lib = CDLL(str(conda_path_unix))
elif conda_path_wind.exists():
    lib = CDLL(str(conda_path_wind))

# fallback: try to load from installed package (site-packages)
if lib is None:
    lib_path = resources.files("pyopentrustregion") / f"libotrtestsuite.{ext}"
    if lib_path.is_file():
        lib = CDLL(str(lib_path))

# fallback: try to load from same directory (editable install)
if lib is None:
    local_path = os.path.join(os.path.dirname(__file__), f"libotrtestsuite.{ext}")
    if os.path.exists(local_path):
        lib = CDLL(local_path)

# fallback: try to load from ../build (development build)
if lib is None:
    build_path = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../build", f"libotrtestsuite.{ext}")
    )
    if os.path.exists(build_path):
        lib = CDLL(build_path)

# if all failed
if lib is None:
    raise FileNotFoundError(
        f"Cannot find any of the expected libraries: libotrtestsuite.{ext}"
    )


# define all tests in alphabetical order
fortran_tests = {
    "opentrustregion_tests": [
        "abs_diag_precond",
        "accept_trust_region_step",
        "add_column",
        "add_error_origin",
        "bisection",
        "bracket",
        "extend_symm_matrix",
        "generate_random_trial_vectors",
        "generate_trial_vectors",
        "gram_schmidt",
        "init_rng",
        "init_solver_settings",
        "init_stability_settings",
        "jacobi_davidson_correction",
        "level_shifted_davidson",
        "level_shifted_diag_precond",
        "print_message",
        "min_eigval",
        "minres",
        "newton_step",
        "orthogonal_projection",
        "print_results",
        "solver",
        "solver_sanity_check",
        "split_string_by_space",
        "stability_check",
        "stability_sanity_check",
        "string_to_lowercase",
        "symm_mat_min_eig",
        "truncated_conjugate_gradient",
    ],
    "c_interface_tests": [
        "assign_solver_c_f",
        "assign_solver_f_c",
        "assign_stability_c_f",
        "assign_stability_f_c",
        "character_from_c",
        "character_to_c",
        "conv_check_f_wrapper",
        "hess_x_f_wrapper",
        "init_solver_settings_c",
        "init_stability_settings_c",
        "logger_f_wrapper",
        "obj_func_f_wrapper",
        "precond_f_wrapper",
        "project_f_wrapper",
        "solver_c_wrapper",
        "stability_check_c_wrapper",
        "update_orbs_f_wrapper",
    ],
    "system_tests": ["h2o_atomic_fb", "h2o_saddle_fb"],
}

# define return type of Fortran functions
for tests in fortran_tests.values():
    for test in tests:
        getattr(lib, f"test_{test}").restype = c_bool


# define function to add tests to test classes
def add_tests(cls):

    def create_test(func_name):
        def test(self):
            result = getattr(lib, func_name)()
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

        # get combined fields
        combined_fields = (
            SolverSettings.c_struct._fields_ + StabilitySettings.c_struct._fields_
        )

        # get non-duplicate fields and allocate memory for reference values
        fields = []
        seen_names = set()
        # loop by type order
        for curr_type in (c_bool, c_real, c_int):
            # solver fields of this type
            for name, t in combined_fields:
                if t == curr_type and name not in seen_names and name != "initialized":
                    fields.append((name + "_ref", t))
                    seen_names.add(name)

        # handle fixed-size c_char arrays (strings)
        for name, t in combined_fields:
            if (
                issubclass(t, Array)
                and getattr(t, "_type_", None) is c_char
                and name not in seen_names
            ):
                fields.append((name + "_ref", t))
                seen_names.add(name)

        # create class to read reference values
        class RefSettingsC(Structure):
            _fields_ = fields

        # create instance
        ref_settings = RefSettingsC()

        # call Fortran to fill values
        lib.get_reference_values.argtypes = [POINTER(RefSettingsC)]
        lib.get_reference_values.restype = None
        lib.get_reference_values(byref(ref_settings))

        # extract Python values
        for name, _ in fields:
            ref_value = getattr(ref_settings, name)
            if isinstance(ref_value, bytes):
                ref_value = ref_value.decode("utf-8")
            setattr(cls, name, ref_value)

        return super().setUpClass()

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.solver", lib.mock_solver)
    def test_solver_py_interface(self):
        """
        this function tests the solver python interface
        """
        n_param = 3

        def mock_obj_func(kappa):
            """
            this function is a mock function for the objective function
            """
            return np.sum(kappa)

        def mock_update_orbs(kappa, grad, h_diag):
            """
            this function is a mock function for the orbital update function
            """
            func = np.sum(kappa)
            grad[:] = 2 * kappa
            h_diag[:] = 3 * kappa

            def hess_x(x, hess_x):
                hess_x[:] = 4 * x

            return func, hess_x

        def mock_precond(residual, mu, precond_residual):
            """
            this function is a mock function for the preconditioner function
            """
            precond_residual[:] = mu * residual

        def mock_project(vector):
            """
            this function is a mock function for the projection function
            """
            vector[:] = 2 * vector

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

        # initialize settings object
        settings = SolverSettings()
        settings.precond = mock_precond
        settings.project = mock_project
        settings.conv_check = mock_conv_check
        settings.logger = mock_logger
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(settings, field_name, getattr(self, field_name + "_ref"))

        # initialize logging boolean
        test_logger = False

        # call solver python interface with optional arguments
        solver(mock_obj_func, mock_update_orbs, n_param, settings)

        # check if logger was called correctly
        if not test_logger:
            print(" test_solver_py_interface failed: Called logging function wrong.")

        self.assertTrue(
            c_bool.in_dll(lib, "test_solver_interface").value and test_logger,
            "test_solver_py_interface failed",
        )
        print(" test_solver_py_interface PASSED")

    # replace original library with mock library
    @patch(
        "pyopentrustregion.python_interface.lib.stability_check",
        lib.mock_stability_check,
    )
    def test_stability_check_py_interface(self):
        """
        this function tests the stability check python interface
        """
        n_param = 3
        h_diag = np.full(n_param, 3.0, dtype=np.float64)

        def mock_hess_x(x, hess_x):
            hess_x[:] = 4 * x

        def mock_precond(residual, mu, precond_residual):
            """
            this function is a mock function for the preconditioner function
            """
            precond_residual[:] = mu * residual

        def mock_project(vector):
            """
            this function is a mock function for the projection function
            """
            vector[:] = 2 * vector

        def mock_logger(message):
            """
            this function is a mock function for the logging function
            """
            nonlocal test_logger
            if message == "test":
                test_logger = True
            return

        # initialize settings object
        settings = StabilitySettings()
        settings.precond = mock_precond
        settings.project = mock_project
        settings.logger = mock_logger
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(settings, field_name, getattr(self, field_name + "_ref"))

        # allocate memory for descent direction
        kappa = np.empty(n_param, dtype=np.float64)

        # initialize logging boolean
        test_logger = False

        # call stability check python interface with optional arguments
        stable = stability_check(h_diag, mock_hess_x, n_param, settings, kappa=kappa)

        # check if logger was called correctly
        if not test_logger:
            print(
                " test_stability_check_py_interface failed: Called logging function "
                "wrong."
            )

        # check if returned variables are correct
        if stable:
            print(
                " test_stability_check_py_interface failed: Returned stability boolean "
                "wrong."
            )
        wrong_direction = not np.allclose(
            kappa, np.full(n_param, 1.0, dtype=np.float64)
        )
        if wrong_direction:
            print(
                " test_stability_check_py_interface failed: Returned direction wrong."
            )

        self.assertTrue(
            c_bool.in_dll(lib, "test_stability_check_interface").value
            and test_logger
            and not stable
            and not wrong_direction,
            "test_stability_check_py_interface failed",
        )
        print(" test_stability_check_py_interface PASSED")

    @patch.object(SolverSettings, "init_c_struct", lib.mock_init_solver_settings)
    def test_solver_settings(self):
        """
        this function ensure the SolverSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = SolverSettings()
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p:
                if (
                    getattr(settings, field_name) is not None
                    or getattr(settings.settings_c, field_name) is not None
                ):
                    print(
                        " test_solver_settings failed: Optional function pointer "
                        f"{field_name} not initialized correctly."
                    )
                    test_passed = False
            elif field_name == "initialized":
                if not getattr(settings, field_name):
                    print(
                        " test_solver_settings failed: Field initialized not "
                        "initialized correctly."
                    )
                    test_passed = False
            else:
                ref_value = getattr(self, field_name + "_ref")
                if field_type == c_real:
                    match = np.isclose(getattr(settings, field_name), ref_value)
                else:
                    match = getattr(settings, field_name) == ref_value
                if not match:
                    print(field_name, getattr(settings, field_name), ref_value)
                    print(
                        f" test_solver_settings failed: Field {field_name} not "
                        "initialized correctly."
                    )
                    test_passed = False

        dummy_error_code = 42

        def dummy_precond():
            return dummy_error_code

        settings.set_optional_callback("precond", CFUNCTYPE(c_int)(dummy_precond))

        c_ptr = getattr(settings.settings_c, "precond")
        c_interface = getattr(settings.settings_c, "precond_interface", None)

        if (
            c_ptr is None
            or (isinstance(c_ptr, c_void_p) and c_ptr.value is None)
            or not callable(c_interface)
            or c_interface() != dummy_error_code
        ):
            print(
                " test_solver_settings failed: Optional callbacks are not set "
                "correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_solver_settings failed")
        print(" test_solver_settings PASSED")

    @patch.object(StabilitySettings, "init_c_struct", lib.mock_init_stability_settings)
    def test_stability_settings(self):
        """
        this function ensure the StabilitySettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = StabilitySettings()
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p:
                if (
                    getattr(settings, field_name) is not None
                    or getattr(settings.settings_c, field_name) is not None
                ):
                    print(
                        " test_stability_settings failed: Optional function pointer "
                        f"{field_name} not initialized correctly."
                    )
                    test_passed = False
            elif field_name == "initialized":
                if not getattr(settings, field_name):
                    print(
                        " test_stability_settings failed: Field initialized not "
                        "initialized correctly."
                    )
                    test_passed = False
            else:
                ref_value = getattr(self, field_name + "_ref")
                if field_type == c_real:
                    match = np.isclose(getattr(settings, field_name), ref_value)
                else:
                    match = getattr(settings, field_name) == ref_value
                if not match:
                    print(field_name, getattr(settings, field_name), ref_value)
                    print(
                        f" test_stability_settings failed: Field {field_name} not "
                        "initialized correctly."
                    )
                    test_passed = False

        self.assertTrue(test_passed, "test_stability_settings failed")
        print(" test_stability_settings PASSED")


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
        test_data = Path(__file__).parent / "test_data"
        if not os.path.isdir(test_data):
            raise RuntimeError(
                "test_data directory does not exist in same directory as testsuite.py."
            )
        lib.set_test_data_path(str(test_data).encode("utf-8"))
        return super().setUpClass()


@unittest.skipUnless(NUMPY_AVAILABLE, "NumPy not available.")
class PySystemTests(unittest.TestCase):
    """
    end-to-end system tests that drive the real pyopentrustregion.solver and
    stability_check wrappers (not the mock library) against the Hartmann 6D
    problem. This is the only test path that exercises the full Python public
    API on a real numerical workload.
    """

    # Hartmann 6D parameters (mirrors tests/opentrustregion_unit_tests.f90)
    alpha = np.array([1.0, 1.2, 3.0, 3.2])
    A = np.array(
        [
            [10.0,  3.0, 17.0,  3.5,  1.7,  8.0],
            [ 0.05, 10.0, 17.0,  0.1,  8.0, 14.0],
            [ 3.0,  3.5,  1.7, 10.0, 17.0,  8.0],
            [17.0,  8.0,  0.05, 10.0,  0.1, 14.0],
        ]
    )
    P = np.array(
        [
            [0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886],
            [0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991],
            [0.2348, 0.1451, 0.3522, 0.2883, 0.3047, 0.6650],
            [0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381],
        ]
    )
    minimum1 = np.array(
        [0.20168951, 0.15001069, 0.47687398,
         0.27533243, 0.31165162, 0.65730053]
    )
    minimum2 = np.array(
        [0.40465313, 0.88244493, 0.84610160,
         0.57398969, 0.13892673, 0.03849589]
    )
    saddle_point = np.array(
        [0.35278250, 0.59374767, 0.47631257,
         0.40058250, 0.31111531, 0.32397158]
    )
    n_param = 6

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running system tests for Python interface...")
        print(50 * "-")
        return super().setUpClass()

    # ---- Hartmann 6D primitives shared by the callbacks ------------------

    @classmethod
    def _exp_terms(cls, x):
        return np.exp(-np.sum(cls.A * (x - cls.P) ** 2, axis=1))

    @classmethod
    def _func(cls, x):
        return -np.dot(cls.alpha, cls._exp_terms(x))

    @classmethod
    def _grad(cls, x):
        e = cls._exp_terms(x)
        return np.array(
            [np.sum(2.0 * cls.alpha * cls.A[:, j] * (x[j] - cls.P[:, j]) * e)
             for j in range(cls.n_param)]
        )

    @classmethod
    def _hess(cls, x):
        e = cls._exp_terms(x)
        H = np.zeros((cls.n_param, cls.n_param))
        for i in range(cls.n_param):
            H[i, i] = 2.0 * np.sum(
                cls.alpha * cls.A[:, i] * e
                * (1.0 - 2.0 * cls.A[:, i] * (x[i] - cls.P[:, i]) ** 2)
            )
            for j in range(i):
                H[i, j] = -4.0 * np.sum(
                    cls.alpha * cls.A[:, i] * cls.A[:, j]
                    * (x[i] - cls.P[:, i]) * (x[j] - cls.P[:, j]) * e
                )
                H[j, i] = H[i, j]
        return H

    # ---- Tests -----------------------------------------------------------

    def test_solver_py_system(self):
        """drive solver() end-to-end on Hartmann 6D from two starting points"""
        # mutable closure state, mirroring the curr_vars module global on the
        # Fortran side
        state = {"curr": np.array([0.20, 0.15, 0.48, 0.28, 0.31, 0.66]),
                 "hess": None,
                 "logger_called": False}

        def update_orbs(delta_vars, grad, h_diag):
            state["curr"] = state["curr"] + delta_vars
            x = state["curr"]
            state["hess"] = self._hess(x)
            grad[:] = self._grad(x)
            h_diag[:] = np.diag(state["hess"])

            def hess_x(v, hv):
                hv[:] = state["hess"] @ v

            return self._func(x), hess_x

        def obj_func(delta_vars):
            return self._func(state["curr"] + delta_vars)

        def precond(residual, mu, out):
            # identity preconditioner; exercises the callback without
            # producing a zero vector when mu=0 (which would trip the
            # Gram-Schmidt zero-vector guard)
            out[:] = residual

        def logger(msg):
            state["logger_called"] = True

        settings = SolverSettings()
        settings.precond = precond
        settings.logger = logger
        settings.verbose = 3  # ensure logger is exercised

        solver(obj_func, update_orbs, self.n_param, settings)
        self.assertTrue(
            np.allclose(state["curr"], self.minimum1, atol=1e-4),
            f"did not converge to minimum1: got {state['curr']}",
        )
        self.assertTrue(state["logger_called"], "logger was never called")

        # restart near the saddle - solver must still reach a known minimum
        state["curr"] = np.array([0.35, 0.59, 0.48, 0.40, 0.31, 0.32])
        solver(obj_func, update_orbs, self.n_param, settings)
        reached_min = (
            np.allclose(state["curr"], self.minimum1, atol=1e-4)
            or np.allclose(state["curr"], self.minimum2, atol=1e-4)
        )
        self.assertTrue(
            reached_min, f"solver from saddle did not reach a minimum: got {state['curr']}"
        )
        print(" test_solver_py_system PASSED")

    def test_stability_check_py_system(self):
        """drive stability_check() end-to-end at a minimum and at a saddle"""
        # at the minimum: must be reported stable
        H = self._hess(self.minimum1)
        h_diag = np.diag(H).copy()

        def hess_x(v, hv):
            hv[:] = H @ v

        settings = StabilitySettings()
        kappa = np.zeros(self.n_param)
        stable = stability_check(h_diag, hess_x, self.n_param, settings, kappa=kappa)
        self.assertTrue(stable, "minimum1 was reported unstable")

        # at the saddle: must be reported unstable, descent direction parallel
        # to the known negative-curvature eigenvector
        H = self._hess(self.saddle_point)
        h_diag = np.diag(H).copy()

        def hess_x_saddle(v, hv):
            hv[:] = H @ v

        kappa = np.zeros(self.n_param)
        stable = stability_check(
            h_diag, hess_x_saddle, self.n_param, settings, kappa=kappa
        )
        self.assertFalse(stable, "saddle_point was reported stable")
        ref = np.array(
            [-0.173375920238, -0.518489821791, -6.432848975252e-3,
             -0.340127852882,  3.066460316955e-3, 0.765095650196]
        )
        self.assertAlmostEqual(
            abs(np.dot(kappa, ref)), 1.0, places=6,
            msg="descent direction does not match the known eigenvector",
        )
        print(" test_stability_check_py_system PASSED")


if __name__ == "__main__":
    unittest.main(verbosity=0)
