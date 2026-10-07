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
    c_char_p,
    c_void_p,
    Array,
    byref,
    CFUNCTYPE,
    Structure,
    POINTER,
)
from unittest.mock import patch
from pathlib import Path
import numpy as np

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
    "opentrustregion_unit_tests": [
        "abs_diag_precond",
        "accept_trust_region_step",
        "add_column",
        "add_error_origin",
        "add_trial_vector",
        "bisection",
        "bracket",
        "check_stationary_point",
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
        "minres",
        "newton_step",
        "option_list",
        "orthogonal_projection",
        "print_message",
        "print_results",
        "solver",
        "solver_sanity_check",
        "split_string_by_space",
        "stability_check",
        "stability_sanity_check",
        "string_to_lowercase",
        "symm_mat_diag",
        "symm_mat_min_eig",
        "truncated_conjugate_gradient",
    ],
    "c_interface_unit_tests": [
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
        "store_optional_c_callbacks",
        "update_orbs_f_wrapper",
    ],
    "opentrustregion_system_tests": [
        "h2o_fb_solver_default",
        "h2o_fb_solver_jacobi_davidson",
        "h2o_fb_solver_line_search",
        "h2o_fb_solver_stability",
        "h2o_fb_solver_tcg",
        "h2o_fb_stability_check_default",
        "h2o_fb_stability_check_jacobi_davidson",
    ],
    "c_integration_tests": [
        "settings_layout",
        "solver_c",
        "solver_settings_init",
        "stability_check_c",
        "stability_settings_init",
    ],
}

# define return type of Fortran functions
for tests in fortran_tests.values():
    for test in tests:
        getattr(lib, f"test_{test}").restype = c_bool

# define signatures of the Fortran functions providing the reference settings
lib.reference_field.restype = None
lib.get_reference_solver_values.argtypes = [
    POINTER(SolverSettings.c_struct),
    c_char_p,
]
lib.get_reference_solver_values.restype = None


def reference(name, field_type):
    """
    this function returns the reference value of a settings field by its name,
    prefixed by "stability_settings." for the nested settings, which also describe
    standalone stability check settings, from the Fortran test reference module
    """
    value = c_real()
    keyword = dict(SolverSettings.c_struct._fields_)["subsystem_solver"]()
    lib.reference_field(name.encode(), byref(value), keyword)
    if issubclass(field_type, Array):
        return keyword.value.decode()
    if field_type == c_bool:
        return bool(value.value)
    if field_type == c_int:
        return int(value.value)
    return value.value


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
class OpenTrustRegionUnitTests(unittest.TestCase):
    """
    this class contains unit tests for opentrustregion
    """

    tests = fortran_tests["opentrustregion_unit_tests"]

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running unit tests for OpenTrustRegion...")
        print(50 * "-")
        return super().setUpClass()


@add_tests
class CInterfaceUnitTests(unittest.TestCase):
    """
    this class contains unit tests for the C interface
    """

    tests = fortran_tests["c_interface_unit_tests"]

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running unit tests for C interface...")
        print(50 * "-")
        return super().setUpClass()


class PyInterfaceUnitTests(unittest.TestCase):
    """
    this class contains unit tests for the Python interface
    """

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running unit tests for Python interface...")
        print(50 * "-")

        return super().setUpClass()

    @staticmethod
    def _mock_precond(residual, mu, precond_residual):
        """
        this function is a mock function for the preconditioner function
        """
        precond_residual[:] = mu * residual

    @staticmethod
    def _mock_project(vector):
        """
        this function is a mock function for the projection function
        """
        vector[:] = 2 * vector

    @staticmethod
    def _raising_logger(message):
        """
        this function is a logging function that raises
        """
        raise ValueError("logging failure")

    # replace original library with mock library
    @patch("pyopentrustregion.python_interface.lib.solver", lib.mock_solver)
    def test_solver_py_interface(self):
        """
        this function tests the solver python interface
        """
        n_param = c_int.in_dll(lib, "test_n_param").value

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

        def mock_conv_check():
            """
            this function is a mock function for the convergence check function
            """
            return True

        # initialize settings object, the logging function records the messages
        messages = []
        settings = SolverSettings()
        settings.precond = self._mock_precond
        settings.project = self._mock_project
        settings.conv_check = mock_conv_check
        settings.logger = messages.append
        settings.stability_settings.precond = self._mock_precond
        settings.stability_settings.project = self._mock_project
        settings.stability_settings.logger = messages.append
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if (
                field_type == c_void_p
                or field_name == "initialized"
                or (isinstance(field_type, type) and issubclass(field_type, Structure))
            ):
                continue
            setattr(settings, field_name, reference(field_name, field_type))

        # set reference values for nested stability check settings
        for field_info in settings.stability_settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(
                settings.stability_settings,
                field_name,
                reference("stability_settings." + field_name, field_type),
            )

        # call solver python interface with optional arguments, the result of the mock
        # is cleared before and read right after the call
        interface_flag = c_bool.in_dll(lib, "test_solver_interface")
        interface_flag.value = False
        solver(mock_obj_func, mock_update_orbs, n_param, settings)
        interface_passed = interface_flag.value

        # check if logger was called correctly
        test_logger = "test" in messages
        if not test_logger:
            print(" test_solver_py_interface failed: Called logging function wrong.")

        # a logging function that raises must not be silent, the exception is reported
        # once the solver has returned
        settings.logger = self._raising_logger
        logger_error_reported = False
        try:
            solver(mock_obj_func, mock_update_orbs, n_param, settings)
        except RuntimeError as e:
            logger_error_reported = isinstance(e.__cause__, ValueError)
        if not logger_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by logging "
                "function was not reported."
            )

        # a callback function that raises makes the solver fail with the exception as
        # the cause of the reported error
        def raising_obj_func(kappa):
            raise ValueError("objective function failure")

        settings.logger = messages.append
        callback_error_reported = False
        try:
            solver(raising_obj_func, mock_update_orbs, n_param, settings)
        except RuntimeError as e:
            callback_error_reported = isinstance(e.__cause__, ValueError)
        if not callback_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by objective "
                "function was not reported."
            )

        self.assertTrue(
            interface_passed
            and test_logger
            and logger_error_reported
            and callback_error_reported,
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
        n_param = c_int.in_dll(lib, "test_n_param").value
        h_diag = np.full(n_param, 3.0, dtype=np.float64)

        def mock_hess_x(x, hess_x):
            hess_x[:] = 4 * x

        # initialize settings object, the logging function records the messages
        messages = []
        settings = StabilitySettings()
        settings.precond = self._mock_precond
        settings.project = self._mock_project
        settings.logger = messages.append
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_type == c_void_p or field_name == "initialized":
                continue
            setattr(
                settings,
                field_name,
                reference("stability_settings." + field_name, field_type),
            )

        # allocate memory for descent direction
        kappa = np.empty(n_param, dtype=np.float64)

        # call stability check python interface with optional arguments, the result of
        # the mock is cleared before and read right after the call
        interface_flag = c_bool.in_dll(lib, "test_stability_check_interface")
        interface_flag.value = False
        stable = stability_check(h_diag, mock_hess_x, n_param, settings, kappa=kappa)
        interface_passed = interface_flag.value

        # check if logger was called correctly
        test_logger = "test" in messages
        if not test_logger:
            print(
                " test_stability_check_py_interface failed: Called logging function "
                "wrong."
            )

        # check if returned variables are correct
        if not stable:
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

        # a logging function that raises must not be silent, the exception is reported
        # once the stability check has returned
        settings.logger = self._raising_logger
        logger_error_reported = False
        try:
            stability_check(h_diag, mock_hess_x, n_param, settings, kappa=kappa)
        except RuntimeError as e:
            logger_error_reported = isinstance(e.__cause__, ValueError)
        if not logger_error_reported:
            print(
                " test_stability_check_py_interface failed: Exception raised by "
                "logging function was not reported."
            )

        # call stability check python interface without a returned direction
        settings.logger = messages.append
        stable_without_direction = stability_check(
            h_diag, mock_hess_x, n_param, settings
        )
        if not stable_without_direction:
            print(
                " test_stability_check_py_interface failed: Returned stability boolean "
                "wrong without returned direction."
            )

        # a callback function that raises makes the stability check fail with the
        # exception as the cause of the reported error
        def raising_hess_x(x, hess_x):
            raise ValueError("Hessian linear transformation failure")

        callback_error_reported = False
        try:
            stability_check(h_diag, raising_hess_x, n_param, settings)
        except RuntimeError as e:
            callback_error_reported = isinstance(e.__cause__, ValueError)
        if not callback_error_reported:
            print(
                " test_stability_check_py_interface failed: Exception raised by "
                "Hessian linear transformation was not reported."
            )

        self.assertTrue(
            interface_passed
            and test_logger
            and stable
            and not wrong_direction
            and logger_error_reported
            and stable_without_direction
            and callback_error_reported,
            "test_stability_check_py_interface failed",
        )
        print(" test_stability_check_py_interface PASSED")

    @staticmethod
    def _settings_hold_reference(test_name, settings, prefix, location=""):
        """
        this function checks that the optional callback functions of a settings object
        are unset and that every other field holds its reference value, the names of
        the stability check settings carry the prefix "stability_settings.", the host
        context is skipped since the Python wrapper never writes it
        """
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_name == "context":
                continue
            if field_type == c_void_p:
                if getattr(settings, field_name) is not None:
                    print(
                        f" {test_name} failed: Optional function pointer {field_name} "
                        f"not initialized correctly{location}."
                    )
                    test_passed = False
            elif isinstance(field_type, type) and issubclass(field_type, Structure):
                test_passed &= PyInterfaceUnitTests._settings_hold_reference(
                    test_name,
                    getattr(settings, field_name),
                    "stability_settings.",
                    " for nested stability settings",
                )
            else:
                ref_value = reference(prefix + field_name, field_type)
                value = getattr(settings, field_name)
                if field_type == c_real:
                    match = np.isclose(value, ref_value)
                else:
                    match = value == ref_value
                if not match:
                    print(
                        f" {test_name} failed: Field {field_name} not initialized "
                        f"correctly{location}."
                    )
                    test_passed = False
        return test_passed

    @patch.object(SolverSettings, "init_c_struct", lib.mock_init_solver_settings)
    def test_solver_settings(self):
        """
        this function ensure the SolverSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = SolverSettings()
        test_passed = self._settings_hold_reference(
            "test_solver_settings", settings, ""
        )

        dummy_error_code = 42

        def dummy_precond():
            return dummy_error_code

        settings.set_optional_callback(
            "precond", dummy_precond, lambda x: x, CFUNCTYPE(c_int)
        )

        c_ptr = getattr(settings.settings_c, "precond")
        c_interface = getattr(settings.settings_c, "precond_interface", None)

        if (
            c_ptr is None
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
        test_passed = self._settings_hold_reference(
            "test_stability_settings", settings, "stability_settings."
        )
        self.assertTrue(test_passed, "test_stability_settings failed")
        print(" test_stability_settings PASSED")


@add_tests
class OpenTrustRegionSystemTests(unittest.TestCase):
    """
    this class contains system tests for opentrustregion
    """

    tests = fortran_tests["opentrustregion_system_tests"]

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


@add_tests
class CIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library through the C
    header
    """

    tests = fortran_tests["c_integration_tests"]

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running integration tests for C interface...")
        print(50 * "-")
        return super().setUpClass()


class PyIntegrationTests(unittest.TestCase):
    """
    this class contains integration tests that drive the real library (not the mock
    library) through the Python interface, with the Hartmann 6D problem as a small
    workload that exercises every callback, settings field and return value crossing
    the interface
    """

    @classmethod
    def setUpClass(cls):
        print(50 * "-")
        print("Running integration tests for Python interface...")
        print(50 * "-")

        # read the Hartmann 6D problem from the library
        def read_array(name, *shape):
            size = int(np.prod(shape))
            array = (c_real * size).in_dll(lib, name)
            return np.frombuffer(array, dtype=np.dtype(c_real), count=size).reshape(
                shape, order="F"
            )

        cls.n_param = c_int.in_dll(lib, "hartmann6d_n_param").value
        cls.n_terms = c_int.in_dll(lib, "hartmann6d_n_terms").value
        cls.alpha = read_array("hartmann6d_alpha", cls.n_terms)
        cls.A = read_array("hartmann6d_A", cls.n_terms, cls.n_param)
        cls.P = read_array("hartmann6d_P", cls.n_terms, cls.n_param)
        cls.minimum1 = read_array("hartmann6d_minimum1", cls.n_param)
        cls.near_minimum = read_array("hartmann6d_near_minimum", cls.n_param)
        cls.saddle_point = read_array("hartmann6d_saddle_point", cls.n_param)

        return super().setUpClass()

    # Hartmann 6D primitives shared by the callbacks

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
            [
                np.sum(2.0 * cls.alpha * cls.A[:, j] * (x[j] - cls.P[:, j]) * e)
                for j in range(cls.n_param)
            ]
        )

    @classmethod
    def _hess(cls, x):
        e = cls._exp_terms(x)
        H = np.zeros((cls.n_param, cls.n_param))
        for i in range(cls.n_param):
            H[i, i] = 2.0 * np.sum(
                cls.alpha
                * cls.A[:, i]
                * e
                * (1.0 - 2.0 * cls.A[:, i] * (x[i] - cls.P[:, i]) ** 2)
            )
            for j in range(i):
                H[i, j] = -4.0 * np.sum(
                    cls.alpha
                    * cls.A[:, i]
                    * cls.A[:, j]
                    * (x[i] - cls.P[:, i])
                    * (x[j] - cls.P[:, j])
                    * e
                )
                H[j, i] = H[i, j]
        return H

    # Tests

    def test_settings_layout(self):
        """
        this function checks that every field of the Python settings structures,
        including the nested stability settings, is read back under its own name from
        settings whose fields Fortran sets one by one to the reference values
        """
        test_passed = True

        # names of the logicals, prefixed for the nested settings
        logicals = [
            name
            for name, field_type in SolverSettings.c_struct._fields_
            if field_type == c_bool
        ] + [
            "stability_settings." + name
            for name, field_type in StabilitySettings.c_struct._fields_
            if field_type == c_bool
        ]

        def read(settings_c, name):
            for part in name.split("."):
                settings_c = getattr(settings_c, part)
            return settings_c

        # check that every logical is read back under its own name, only one is set at
        # a time so that swapped logicals can be told apart
        for true_logical in logicals:
            settings_c = SolverSettings.c_struct()
            lib.get_reference_solver_values(byref(settings_c), true_logical.encode())
            for name in logicals:
                if read(settings_c, name) != (name == true_logical):
                    print(f" test_settings_layout failed: Field {name} misplaced.")
                    test_passed = False

        # check that every other field is read back under its own name
        settings_c = SolverSettings.c_struct()
        lib.get_reference_solver_values(byref(settings_c), None)
        for struct, prefix in [
            (settings_c, ""),
            (settings_c.stability_settings, "stability_settings."),
        ]:
            for name, field_type in struct._fields_:
                if field_type == c_bool or field_type == StabilitySettings.c_struct:
                    continue
                value = getattr(struct, name)
                if field_type == c_void_p:
                    value = value or 0
                elif issubclass(field_type, Array):
                    value = value.decode()
                if value != reference(prefix + name, field_type):
                    print(
                        f" test_settings_layout failed: Field {prefix + name} "
                        "misplaced."
                    )
                    test_passed = False
        self.assertTrue(test_passed, "test_settings_layout failed")
        print(" test_settings_layout PASSED")

    def test_solver_py(self):
        """
        this function drives solver() end-to-end from near a minimum with every
        optional callback set
        """
        test_passed = True

        # mutable closure state holding the current point, starting in the quadratic
        # region near the first minimum, its Hessian and the names of the callbacks
        # that were reached
        state = {
            "curr": self.near_minimum.copy(),
            "hess": None,
            "called": set(),
        }

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

        # identity preconditioners and projections and loggers that record their
        # names, the stability variants are set only on the nested stability settings,
        # so the test can tell whether the internal stability check reached its own
        # callback slots rather than the solver's, the identity preconditioners
        # exercise the callback without producing a zero vector when mu=0 (which would
        # trip the Gram-Schmidt zero-vector guard)
        def recording_precond(name):
            def precond(residual, mu, out):
                state["called"].add(name)
                out[:] = residual

            return precond

        def recording_callback(name):
            def callback(*args):
                state["called"].add(name)

            return callback

        def conv_check():
            # never report convergence, so the solve is unaffected
            state["called"].add("conv_check")
            return False

        settings = SolverSettings()
        settings.precond = recording_precond("precond")
        settings.project = recording_callback("project")
        settings.conv_check = conv_check
        settings.logger = recording_callback("logger")
        settings.stability = True
        settings.stability_settings.precond = recording_precond("stability_precond")
        settings.stability_settings.project = recording_callback("stability_project")
        settings.stability_settings.logger = recording_callback("stability_logger")
        settings.verbose = 3  # ensure the loggers are exercised

        # the maximum precision flag starts opposite to the one the converging solve
        # returns, so that its write-back is detected
        settings.max_precision_reached = True

        solver(obj_func, update_orbs, self.n_param, settings)
        if settings.max_precision_reached:
            print(
                " test_solver_py failed: Maximum precision reached flag was not "
                "returned."
            )
            test_passed = False
        for name, description in [
            ("precond", "Preconditioner"),
            ("project", "Projection"),
            ("conv_check", "Convergence check"),
            ("logger", "Logger"),
        ]:
            if name not in state["called"]:
                print(f" test_solver_py failed: {description} was not called.")
                test_passed = False
        if settings.n_update_orbs <= 0:
            print(" test_solver_py failed: Orbital update counter was not populated.")
            test_passed = False
        if settings.n_hess_x <= 0:
            print(
                " test_solver_py failed: Hessian linear transformation counter was not "
                "populated."
            )
            test_passed = False
        if settings.stability_settings.n_hess_x <= 0:
            print(
                " test_solver_py failed: Hessian linear transformation counter of the "
                "internal stability check was not populated."
            )
            test_passed = False
        for name, description in [
            ("stability_precond", "Preconditioner"),
            ("stability_project", "Projection"),
            ("stability_logger", "Logger"),
        ]:
            if name not in state["called"]:
                print(
                    f" test_solver_py failed: {description} set on the nested "
                    "stability settings was not called by the internal stability check."
                )
                test_passed = False
        self.assertTrue(test_passed, "test_solver_py failed")
        print(" test_solver_py PASSED")

    def test_stability_check_py(self):
        """
        this function drives stability_check() end-to-end at a minimum and at a saddle
        point with every optional callback set
        """
        test_passed = True
        called = set()

        # identity preconditioner, identity projection and logger that record their
        # names, the identity preconditioner exercises the callback without producing
        # a zero vector when mu=0 (which would trip the Gram-Schmidt zero-vector guard)
        def precond(residual, mu, out):
            called.add("precond")
            out[:] = residual

        def project(vector):
            called.add("project")

        def logger(msg):
            called.add("logger")

        settings = StabilitySettings()
        settings.precond = precond
        settings.project = project
        settings.logger = logger
        settings.verbose = 3  # ensure logger is exercised

        # at a minimum, expect stable
        H = self._hess(self.minimum1)
        h_diag = np.diag(H).copy()

        def hess_x(v, hv):
            hv[:] = H @ v

        kappa = np.zeros(self.n_param)
        stable = stability_check(h_diag, hess_x, self.n_param, settings, kappa=kappa)
        if not stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of minimum."
            )
            test_passed = False
        for name, description in [
            ("precond", "Preconditioner"),
            ("project", "Projection"),
            ("logger", "Logger"),
        ]:
            if name not in called:
                print(f" test_stability_check_py failed: {description} was not called.")
                test_passed = False
        if settings.n_hess_x <= 0:
            print(
                " test_stability_check_py failed: Hessian linear transformation "
                "counter was not populated."
            )
            test_passed = False

        # also exercise the call without a direction at the minimum, where the
        # stability flag differs from the False the wrapper starts from
        stable = stability_check(h_diag, hess_x, self.n_param, settings)
        if not stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of minimum when not passing direction."
            )
            test_passed = False

        # at a saddle, the returned direction replaces the zero direction returned at
        # the minimum and has to be a normalized direction of negative curvature
        H = self._hess(self.saddle_point)
        h_diag = np.diag(H).copy()
        stability_check(h_diag, hess_x, self.n_param, settings, kappa=kappa)
        if abs(np.linalg.norm(kappa) - 1.0) > 1e-6:
            print(
                " test_stability_check_py failed: Stability check does not return a "
                "normalized direction for saddle point."
            )
            test_passed = False
        if kappa @ H @ kappa >= 0.0:
            print(
                " test_stability_check_py failed: Stability check does not return a "
                "direction of negative curvature for saddle point."
            )
            test_passed = False
        self.assertTrue(test_passed, "test_stability_check_py failed")
        print(" test_stability_check_py PASSED")


if __name__ == "__main__":
    unittest.main(verbosity=0)
