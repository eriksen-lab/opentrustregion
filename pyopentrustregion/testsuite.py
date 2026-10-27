# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import sys
import copy
import signal
import threading
import unittest
import warnings
from ctypes import (
    CDLL,
    c_bool,
    pointer,
    c_char_p,
    c_void_p,
    sizeof,
    memset,
    addressof,
    Array,
    byref,
    cast,
    CFUNCTYPE,
    Structure,
    POINTER,
)
from unittest.mock import patch
from functools import partial
from types import SimpleNamespace
from pathlib import Path
import numpy as np

# when run as a script, import the package this file belongs to rather than another
# one that is installed
if not __package__:
    sys.path.insert(0, str(Path(__file__).parent.parent.absolute()))
from pyopentrustregion.python_interface import (
    SolverSettings,
    StabilitySettings,
    SolverState,
    StabilityState,
    solver,
    stability_check,
    c_int,
    c_real,
    ext,
    _load_lib,
    interruption_state,
    defer_interruption,
    ObjFuncInterface,
    HessXInterface,
    UpdateOrbsInterface,
    PrecondInterface,
    ProjectInterface,
    GetExtraTrialVectorsInterface,
    ConvCheckInterface,
    LoggerInterface,
    raise_on_failure,
    call_library,
    check_struct_sizes,
    hess_x_interface_type,
    precond_interface_type,
    project_interface_type,
    get_extra_trial_vectors_interface_type,
    conv_check_interface_type,
)

# load the testsuite library
lib = _load_lib([f"libotrtestsuite.{ext}"])


# define all tests in alphabetical order
fortran_tests = {
    "opentrustregion_unit_tests": [
        "abs_diag_precond",
        "accept_trust_region_step",
        "add_column",
        "add_error_origin",
        "add_trial_vector",
        "bisection_ah",
        "bisection_mu",
        "block_davidson",
        "bracket",
        "check_stationary_point",
        "extend_matrix",
        "general_mat_diag",
        "generate_random_trial_vectors",
        "generate_trial_vectors",
        "get_stability_trial_space",
        "gram_schmidt",
        "init_defaults",
        "init_rng",
        "jacobi_davidson_correction",
        "level_shifted_davidson",
        "level_shifted_diag_precond",
        "lowest_h_diag_unit_vectors",
        "mat_diag",
        "mat_lowest_eigpairs",
        "mat_min_eig",
        "micro_conv_tol",
        "minres",
        "newton_step",
        "noise_level",
        "option_list",
        "orthogonal_projection",
        "orthogonalize_trial_vectors",
        "perturb_vector",
        "print_message",
        "print_results",
        "saddle_point_step",
        "solver",
        "solver_sanity_check",
        "split_string_by_space",
        "stability_check",
        "stability_sanity_check",
        "string_in",
        "string_to_lowercase",
        "symm_mat_diag",
        "truncated_conjugate_gradient",
    ],
    "c_interface_unit_tests": [
        "assign_solver_c_f",
        "assign_solver_f_c",
        "assign_solver_state_c_f",
        "assign_stability_c_f",
        "assign_stability_f_c",
        "assign_stability_state_c_f",
        "character_from_c",
        "character_to_c",
        "conv_check_f_wrapper",
        "copy_state_to_c",
        "get_extra_trial_vectors_f_wrapper",
        "hess_x_f_wrapper",
        "init_solver_settings_c",
        "init_stability_settings_c",
        "logger_f_wrapper",
        "obj_func_f_wrapper",
        "precond_f_wrapper",
        "project_f_wrapper",
        "solver_c_wrapper",
        "solver_settings_size",
        "solver_state_size",
        "stability_check_c_wrapper",
        "stability_settings_size",
        "stability_state_size",
        "store_optional_c_callbacks",
        "update_orbs_f_wrapper",
    ],
    "opentrustregion_system_tests": [
        "h2o_fb_solver_davidson_ah",
        "h2o_fb_solver_default",
        "h2o_fb_solver_hess_symm",
        "h2o_fb_solver_jacobi_davidson_ah",
        "h2o_fb_solver_jacobi_davidson_ls",
        "h2o_fb_solver_line_search",
        "h2o_fb_solver_stability",
        "h2o_fb_solver_tcg",
        "h2o_fb_stability_check_default",
        "h2o_fb_stability_check_hess_symm",
        "h2o_fb_stability_check_jacobi_davidson",
    ],
    "c_integration_tests": [
        "settings_layout",
        "solver_c",
        "solver_settings_init",
        "stability_check_c",
        "stability_settings_init",
        "state_layout",
    ],
}

# define signatures of the Fortran functions providing the reference settings
lib.reference_field.restype = None
lib.get_reference_solver_values.argtypes = [
    POINTER(SolverSettings.c_struct),
    c_char_p,
]
lib.get_reference_solver_values.restype = None

# define signatures of the Fortran functions providing the reference states
lib.reference_state_field.argtypes = [c_char_p]
lib.reference_state_field.restype = c_real
lib.get_reference_solver_state_values.argtypes = [POINTER(SolverState.c_struct)]
lib.get_reference_solver_state_values.restype = None
lib.get_reference_stability_state_values.argtypes = [POINTER(StabilityState.c_struct)]
lib.get_reference_stability_state_values.restype = None


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


def reference_state(name, field_type):
    """
    this function returns the reference value of a state field by its name, prefixed
    by "stability_state." for the nested state, which also describes a standalone
    stability check state, from the Fortran test reference module
    """
    value = lib.reference_state_field(name.encode())
    if field_type == c_bool:
        return bool(value)
    if field_type == c_int:
        return int(value)
    return value


# define function to add tests to test classes
def add_tests(cls):

    def create_test(func_name):
        def test(self):
            func = getattr(lib, func_name)
            func.restype = c_bool
            result = func()
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
    def _mock_get_extra_trial_vectors(trial_vectors):
        """
        this function is a mock function for the extra trial vector function, which
        fills the i-th trial vector with i
        """
        trial_vectors[:] = np.arange(1, len(trial_vectors) + 1).reshape(-1, 1)

    @staticmethod
    def _mock_conv_check():
        """
        this function is a mock function for the convergence check function
        """
        return True

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

        # the callback functions record that they were called and the logging
        # functions record the messages, separately for the settings and the nested
        # settings so that their callback functions can be told apart
        called = set()

        def recording(name, func):
            def callback(*args):
                called.add(name)
                return func(*args)

            return callback

        messages = []
        stability_messages = []
        settings = SolverSettings()
        settings.precond = recording("precond", self._mock_precond)
        settings.project = recording("project", self._mock_project)
        settings.get_extra_trial_vectors = recording(
            "get_extra_trial_vectors", self._mock_get_extra_trial_vectors
        )
        settings.conv_check = self._mock_conv_check
        settings.logger = messages.append
        settings.stability_settings.precond = recording(
            "stability_precond", self._mock_precond
        )
        settings.stability_settings.project = recording(
            "stability_project", self._mock_project
        )
        settings.stability_settings.get_extra_trial_vectors = recording(
            "stability_get_extra_trial_vectors", self._mock_get_extra_trial_vectors
        )
        settings.stability_settings.conv_check = recording(
            "stability_conv_check", self._mock_conv_check
        )
        settings.stability_settings.logger = stability_messages.append
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
        # is cleared before and read right after the call, and check that the call
        # leaves the settings object unchanged, so that calls running at the same time
        # can share it
        interface_flag = c_bool.in_dll(lib, "test_solver_interface")
        interface_flag.value = False
        state = SolverState()
        settings_c_before = bytes(settings.settings_c)
        solver(mock_update_orbs, mock_obj_func, n_param, settings, state=state)
        interface_passed = interface_flag.value
        settings_kept = bytes(settings.settings_c) == settings_c_before
        if not settings_kept:
            print(
                " test_solver_py_interface failed: Settings object written by the call."
            )

        # check that the state is passed, which the mock fills with the reference
        # values
        state_passed = self._state_holds_reference(
            "test_solver_py_interface", state, ""
        )

        # check if every callback function was called through its own slot
        test_logger = "test" in messages and "test" in stability_messages
        if not test_logger:
            print(" test_solver_py_interface failed: Called logging function wrong.")
        callbacks_called = True
        for name in [
            "precond",
            "project",
            "get_extra_trial_vectors",
            "stability_precond",
            "stability_project",
            "stability_get_extra_trial_vectors",
            "stability_conv_check",
        ]:
            if name not in called:
                print(
                    f" test_solver_py_interface failed: Callback function {name} was "
                    "not called."
                )
                callbacks_called = False

        # a logging function that raises must not be silent, the exception is reported
        # as the cause of a warning once the solver has returned without an error
        settings.logger = self._raising_logger
        logger_error_reported = False
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                solver(mock_update_orbs, mock_obj_func, n_param, settings)
                logger_error_reported = any(
                    issubclass(w.category, RuntimeWarning)
                    and isinstance(w.message.__cause__, ValueError)
                    for w in caught
                )
            except RuntimeError:
                pass
        if not logger_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by logging "
                "function was not reported as a warning."
            )

        # a callback function that raises makes the solver fail with the exception as
        # the cause of the reported error
        def raising_obj_func(kappa):
            raise ValueError("objective function failure")

        settings.logger = messages.append
        callback_error_reported = False
        try:
            solver(mock_update_orbs, raising_obj_func, n_param, settings)
        except RuntimeError as e:
            callback_error_reported = isinstance(e.__cause__, ValueError)
        if not callback_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by objective "
                "function was not reported."
            )

        # an orbital update that raises makes the solver fail with the exception as the
        # cause of the reported error, which a library calling it right away shows
        update_orbs_failure = ValueError("orbital update failure")

        def raising_update_orbs(kappa, grad, h_diag):
            raise update_orbs_failure

        def library_calling_update_orbs(update_orbs_interface, *args):
            vector = np.zeros(n_param)
            vector_ptr = vector.ctypes.data_as(POINTER(c_real))
            return update_orbs_interface(
                vector_ptr,
                pointer(c_real()),
                vector_ptr,
                vector_ptr,
                pointer(hess_x_interface_type()),
                None,
            )

        update_orbs_error_reported = False
        with patch(
            "pyopentrustregion.python_interface.lib.solver",
            library_calling_update_orbs,
        ):
            try:
                solver(raising_update_orbs, mock_obj_func, n_param)
            except RuntimeError as e:
                update_orbs_error_reported = e.__cause__ is update_orbs_failure
        if not update_orbs_error_reported:
            print(
                " test_solver_py_interface failed: Exception raised by orbital update "
                "was not reported."
            )

        # the solver marks the library as innermost while it runs, around the callback
        # function marked within it, and raises an interruption of a callback function
        # as it is
        interruption = KeyboardInterrupt()
        marks = []

        def interrupted_obj_func(kappa):
            marks.append(interruption_state.stack)
            raise interruption

        raised = None
        try:
            solver(mock_update_orbs, interrupted_obj_func, n_param, settings)
        except KeyboardInterrupt as e:
            raised = e
        interruption_raised = raised is interruption and marks == [(False, True)]
        if not interruption_raised:
            print(
                " test_solver_py_interface failed: Library not marked as innermost or "
                "interruption of callback function not raised."
            )

        # call solver python interface without settings, which have to be passed as
        # none so that the library uses its default settings, the results of the mock
        # start opposite to the expected ones
        settings_flag = c_bool.in_dll(lib, "solver_settings_passed")
        settings_flag.value = True
        interface_flag.value = False
        solver(mock_update_orbs, mock_obj_func, n_param)
        no_settings_passed = interface_flag.value and not settings_flag.value
        if settings_flag.value:
            print(
                " test_solver_py_interface failed: Settings passed although none were "
                "given."
            )

        self.assertTrue(
            interface_passed
            and settings_kept
            and state_passed
            and test_logger
            and callbacks_called
            and logger_error_reported
            and callback_error_reported
            and update_orbs_error_reported
            and interruption_raised
            and no_settings_passed,
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
        settings.get_extra_trial_vectors = self._mock_get_extra_trial_vectors
        settings.conv_check = self._mock_conv_check
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
        kappa = np.zeros(n_param, dtype=np.float64)

        # call stability check python interface with optional arguments, the result of
        # the mock is cleared before and read right after the call, and check that the
        # call leaves the settings object unchanged, so that calls running at the same
        # time can share it
        interface_flag = c_bool.in_dll(lib, "test_stability_check_interface")
        interface_flag.value = False
        state = StabilityState()
        settings_c_before = bytes(settings.settings_c)
        stable = stability_check(
            h_diag, mock_hess_x, n_param, settings, kappa=kappa, state=state
        )
        interface_passed = interface_flag.value
        settings_kept = bytes(settings.settings_c) == settings_c_before
        if not settings_kept:
            print(
                " test_stability_check_py_interface failed: Settings object written by "
                "the call."
            )

        # check that the state is passed, which the mock fills with the reference
        # values
        state_passed = self._state_holds_reference(
            "test_stability_check_py_interface", state, "stability_state."
        )

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
        # as the cause of a warning once the stability check has returned without an
        # error, which still returns its result
        settings.logger = self._raising_logger
        logger_error_reported = False
        stable_with_failed_logger = False
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                stable_with_failed_logger = stability_check(
                    h_diag, mock_hess_x, n_param, settings, kappa=kappa
                )
                logger_error_reported = any(
                    issubclass(w.category, RuntimeWarning)
                    and isinstance(w.message.__cause__, ValueError)
                    for w in caught
                )
            except RuntimeError:
                pass
        if not logger_error_reported:
            print(
                " test_stability_check_py_interface failed: Exception raised by "
                "logging function was not reported as a warning."
            )
        if not stable_with_failed_logger:
            print(
                " test_stability_check_py_interface failed: Returned stability boolean "
                "wrong with failed logging function."
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

        # the stability check marks the library as innermost while it runs, around the
        # callback function marked within it, and raises an interruption of a callback
        # function as it is
        interruption = KeyboardInterrupt()
        marks = []

        def interrupted_hess_x(x, hess_x):
            marks.append(interruption_state.stack)
            raise interruption

        raised = None
        try:
            stability_check(h_diag, interrupted_hess_x, n_param, settings)
        except KeyboardInterrupt as e:
            raised = e
        interruption_raised = raised is interruption and marks == [(False, True)]
        if not interruption_raised:
            print(
                " test_stability_check_py_interface failed: Library not marked as "
                "innermost or interruption of callback function not raised."
            )

        # a Hessian diagonal given as a strided view, of another type or as a list is
        # passed with its values, while a Hessian diagonal of another length and an
        # array for the direction which the library can not write into are refused
        matrix = np.zeros((n_param, n_param))
        np.fill_diagonal(matrix, 3.0)
        h_diag_converted = True
        for h_diag_given, case in [
            (np.diag(matrix), "a strided view"),
            (np.full(n_param, 3.0, dtype=np.float32), "a single precision array"),
            ([3.0] * n_param, "a list"),
        ]:
            interface_flag.value = False
            stability_check(h_diag_given, mock_hess_x, n_param, settings)
            if not interface_flag.value:
                print(
                    " test_stability_check_py_interface failed: Hessian diagonal "
                    f"given as {case} not passed correctly."
                )
                h_diag_converted = False
        read_only = np.zeros(n_param)
        read_only.flags.writeable = False
        wrong_inputs_refused = True
        for h_diag_given, kappa_given, case in [
            (np.full(n_param + 1, 3.0), None, "Hessian diagonal of another length"),
            (h_diag, np.zeros(n_param + 1), "Direction of another length"),
            (h_diag, np.zeros(n_param, dtype=np.float32), "Single precision direction"),
            (h_diag, np.zeros(2 * n_param)[::2], "Strided direction"),
            (h_diag, read_only, "Read-only direction"),
            (h_diag, [0.0] * n_param, "Direction given as a list"),
        ]:
            try:
                stability_check(
                    h_diag_given, mock_hess_x, n_param, settings, kappa=kappa_given
                )
                print(f" test_stability_check_py_interface failed: {case} not refused.")
                wrong_inputs_refused = False
            except ValueError:
                pass

        # call stability check python interface without settings, which have to be
        # passed as none so that the library uses its default settings, the results of
        # the mock start opposite to the expected ones
        settings_flag = c_bool.in_dll(lib, "stability_settings_passed")
        settings_flag.value = True
        interface_flag.value = False
        stable_without_settings = stability_check(h_diag, mock_hess_x, n_param)
        no_settings_passed = (
            interface_flag.value and not settings_flag.value and stable_without_settings
        )
        if settings_flag.value:
            print(
                " test_stability_check_py_interface failed: Settings passed although "
                "none were given."
            )

        self.assertTrue(
            interface_passed
            and settings_kept
            and state_passed
            and test_logger
            and stable
            and not wrong_direction
            and logger_error_reported
            and stable_with_failed_logger
            and stable_without_direction
            and callback_error_reported
            and interruption_raised
            and h_diag_converted
            and wrong_inputs_refused
            and no_settings_passed,
            "test_stability_check_py_interface failed",
        )
        print(" test_stability_check_py_interface PASSED")

    @staticmethod
    def _settings_hold_reference(test_name, settings, prefix, location=""):
        """
        this function checks that the optional callback functions of a settings object
        are unset and that every other field holds its reference value, the names of
        the stability check settings carry the prefix "stability_settings.", the host
        context and the initialization flag are skipped since the Python wrapper
        exposes neither
        """
        test_passed = True
        for field_info in settings.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if field_name in ("context", "initialized"):
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

    def _settings_safeguarded(self, test_name, settings):
        """
        this function checks that a settings object accepts no attribute that is not a
        setting or a callback function and that its copies hold the settings in a C
        struct of their own, a copy sharing the callback functions and a deep copy
        copying them, also for the nested stability check settings
        """
        test_passed = True

        # neither the initialization flag nor the host context nor a misspelled
        # setting can be set
        for name in ["initialized", "context", "conv_tl"]:
            try:
                setattr(settings, name, 0)
                print(f" {test_name} failed: Attribute {name} could be set.")
                test_passed = False
            except AttributeError:
                pass

        # change a setting, so that only a copy of the values holds it, and set a
        # callback function that a deep copy copies, on the nested settings too
        all_settings = [settings]
        if hasattr(settings, "stability_settings"):
            all_settings.append(settings.stability_settings)
        for some_settings in all_settings:
            some_settings.conv_tol = 3 * some_settings.conv_tol
            some_settings.precond = partial(self._mock_precond)

        for copy_function, shared, case in [
            (copy.copy, True, "copy"),
            (copy.deepcopy, False, "deep copy"),
        ]:
            settings_copy = copy_function(settings)
            copies = [settings_copy]
            if hasattr(settings_copy, "stability_settings"):
                copies.append(settings_copy.stability_settings)
            for original, duplicate in zip(all_settings, copies):
                # the copy holds the same settings, a change of which leaves the
                # original unchanged
                if bytes(duplicate.settings_c) != bytes(original.settings_c):
                    print(f" {test_name} failed: The {case} holds other settings.")
                    test_passed = False
                conv_tol = original.conv_tol
                duplicate.conv_tol = 2 * conv_tol
                if original.conv_tol != conv_tol:
                    print(
                        f" {test_name} failed: Changing the {case} changed the "
                        "original."
                    )
                    test_passed = False

                # a copy shares the callback functions, a deep copy copies them
                if (duplicate.precond is original.precond) != shared or not callable(
                    duplicate.precond
                ):
                    print(
                        f" {test_name} failed: Callback function wrong for the "
                        f"{case}."
                    )
                    test_passed = False

        return test_passed

    @staticmethod
    def _callbacks_wired(test_name, settings_class):
        """
        this function checks that the interfaces a settings object sets up for a call
        report an exception raised by an optional callback function as an error and
        keep it with the exceptions of the call, also for the nested stability check
        settings
        """
        test_passed = True
        n_param = 2
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))
        callbacks = {
            "precond": (
                precond_interface_type,
                (vector_ptr, pointer(c_real(1.0)), vector_ptr),
            ),
            "project": (project_interface_type, (vector_ptr,)),
            "get_extra_trial_vectors": (
                get_extra_trial_vectors_interface_type,
                (vector_ptr, 1),
            ),
            "conv_check": (conv_check_interface_type, (pointer(c_bool()),)),
        }

        def raising(failure):
            def callback(*args):
                raise failure

            return callback

        # every optional callback function raises an exception of its own
        settings = settings_class()
        all_settings = [(settings, "")]
        if hasattr(settings, "stability_settings"):
            all_settings.append((settings.stability_settings, "stability_settings."))
        failures = {}
        for some_settings, prefix in all_settings:
            for name in callbacks:
                failures[prefix + name] = ValueError(prefix + name)
                setattr(some_settings, name, raising(failures[prefix + name]))

        # set up the interfaces in a copy of the C struct and call each of them
        settings_c = settings.c_struct.from_buffer_copy(settings.settings_c)
        exception = {}
        interfaces = []
        settings.set_optional_callbacks(settings_c, interfaces, n_param, exception)
        all_settings_c = [(settings_c, "")]
        if hasattr(settings, "stability_settings"):
            all_settings_c.append(
                (settings_c.stability_settings, "stability_settings.")
            )
        for some_settings_c, prefix in all_settings_c:
            for name, (interface_type, args) in callbacks.items():
                exception.clear()
                interface = cast(
                    c_void_p(getattr(some_settings_c, name)), interface_type
                )
                if interface(*args, None) != 1 or (
                    exception.get("exc") is not failures[prefix + name]
                ):
                    print(
                        f" {test_name} failed: Exception raised by {prefix}{name} not "
                        "reported as error or not kept."
                    )
                    test_passed = False

        return test_passed

    @staticmethod
    def _state_holds_reference(test_name, state, prefix, location=""):
        """
        this function checks that every field of a state object holds its reference
        value, the names of the stability check state carry the prefix
        "stability_state."
        """
        test_passed = True
        for field_info in state.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if isinstance(field_type, type) and issubclass(field_type, Structure):
                test_passed &= PyInterfaceUnitTests._state_holds_reference(
                    test_name,
                    getattr(state, field_name),
                    "stability_state.",
                    " for nested stability check state",
                )
            elif getattr(state, field_name) != reference_state(
                prefix + field_name, field_type
            ):
                print(
                    f" {test_name} failed: Field {field_name} not read correctly"
                    f"{location}."
                )
                test_passed = False
        return test_passed

    def _state_copies_kept(self, test_name, state, prefix):
        """
        this function checks that a copy and a deep copy of a state holding the
        reference values keep them when the original state is reset, since they hold
        them in a C struct of their own
        """
        test_passed = True
        copies = [copy.copy(state), copy.deepcopy(state)]
        memset(addressof(state.state_c), 0, sizeof(state.c_struct))
        for state_copy, case in zip(copies, ["copy", "deep copy"]):
            test_passed &= self._state_holds_reference(
                test_name, state_copy, prefix, f" for the {case} of a reset state"
            )
        return test_passed

    @staticmethod
    def _state_holds_defaults(test_name, state, location=""):
        """
        this function checks that every field of a state object holds its default
        value, zero or false
        """
        test_passed = True
        for field_info in state.c_struct._fields_:
            field_name, field_type = field_info[:2]
            if isinstance(field_type, type) and issubclass(field_type, Structure):
                test_passed &= PyInterfaceUnitTests._state_holds_defaults(
                    test_name,
                    getattr(state, field_name),
                    " for nested stability check state",
                )
            elif getattr(state, field_name) != 0:
                print(
                    f" {test_name} failed: Field {field_name} does not start from "
                    f"its default value{location}."
                )
                test_passed = False
        return test_passed

    def test_defer_interruption(self):
        """
        this function tests that the handler of SIGINT only records an interruption
        while the library runs and raises it otherwise
        """
        test_passed = True
        try:
            for stack, deferred, case in [
                ((False,), True, "while the library is innermost"),
                ((False, True, False), True, "while a nested call is innermost"),
                ((False, True), False, "while a callback function is innermost"),
                ((), False, "outside a call into the library"),
            ]:
                interruption_state.stack = stack
                interruption_state.pending = None
                raised = False
                try:
                    defer_interruption(signal.SIGINT, None)
                except KeyboardInterrupt:
                    raised = True
                if raised == deferred:
                    print(
                        " test_defer_interruption failed: Interruption "
                        f"{'raised' if raised else 'not raised'} {case}."
                    )
                    test_passed = False
                if (
                    isinstance(interruption_state.pending, KeyboardInterrupt)
                    != deferred
                ):
                    print(
                        " test_defer_interruption failed: Interruption "
                        f"{'not recorded' if deferred else 'recorded'} {case}."
                    )
                    test_passed = False
        finally:
            interruption_state.stack = ()
            interruption_state.pending = None

        self.assertTrue(test_passed, "test_defer_interruption failed")
        print(" test_defer_interruption PASSED")

    def test_callback_interfaces(self):
        """
        this function tests that the interfaces to the Python callback functions
        report an exception or an interruption as an error and keep it, also one of the
        interface itself, that they mark a callback function as innermost while it
        runs, that the convergence check interface passes on the result and that the
        logging interface keeps its first exception and records an interruption
        """
        test_passed = True
        n_param = 2
        vector = np.zeros(n_param)
        vector_ptr = vector.ctypes.data_as(POINTER(c_real))

        def raising(error_type):
            def callback(*args):
                raise error_type("callback function failure")

            return callback

        for error_type in (ValueError, KeyboardInterrupt):
            for name, make_interface, args in [
                (
                    "objective function",
                    lambda exc: ObjFuncInterface(raising(error_type), n_param, exc),
                    (vector_ptr, pointer(c_real())),
                ),
                (
                    "Hessian linear transformation",
                    lambda exc: HessXInterface(raising(error_type), n_param, exc),
                    (vector_ptr, vector_ptr),
                ),
                (
                    "orbital update",
                    lambda exc: UpdateOrbsInterface(raising(error_type), n_param, exc),
                    (vector_ptr, pointer(c_real()), vector_ptr, vector_ptr, None),
                ),
                (
                    "preconditioner",
                    lambda exc: PrecondInterface(raising(error_type), n_param, exc),
                    (vector_ptr, pointer(c_real(1.0)), vector_ptr),
                ),
                (
                    "projection",
                    lambda exc: ProjectInterface(raising(error_type), n_param, exc),
                    (vector_ptr,),
                ),
                (
                    "extra trial vector function",
                    lambda exc: GetExtraTrialVectorsInterface(
                        raising(error_type), n_param, exc
                    ),
                    (vector_ptr, 1),
                ),
                (
                    "convergence check",
                    lambda exc: ConvCheckInterface(raising(error_type), exc),
                    (pointer(c_bool()),),
                ),
            ]:
                exception = {}
                if make_interface(exception)(*args, None) != 1:
                    print(
                        f" test_callback_interfaces failed: {error_type.__name__} "
                        f"raised by {name} not reported as error."
                    )
                    test_passed = False
                if not isinstance(exception.get("exc"), error_type):
                    print(
                        f" test_callback_interfaces failed: {error_type.__name__} "
                        f"raised by {name} not kept."
                    )
                    test_passed = False

        # the interfaces hand the callback functions their inputs read-only, so that
        # writing into one is reported as an error of the callback function and leaves
        # the input unchanged
        def writing_into_input(array, *args):
            array[0] = 1.0

        output = np.zeros(n_param)
        output_ptr = output.ctypes.data_as(POINTER(c_real))
        for name, make_interface, args in [
            (
                "objective function",
                lambda exc: ObjFuncInterface(writing_into_input, n_param, exc),
                (vector_ptr, pointer(c_real())),
            ),
            (
                "Hessian linear transformation",
                lambda exc: HessXInterface(writing_into_input, n_param, exc),
                (vector_ptr, output_ptr),
            ),
            (
                "orbital update",
                lambda exc: UpdateOrbsInterface(writing_into_input, n_param, exc),
                (vector_ptr, pointer(c_real()), output_ptr, output_ptr, None),
            ),
            (
                "preconditioner",
                lambda exc: PrecondInterface(writing_into_input, n_param, exc),
                (vector_ptr, pointer(c_real(1.0)), output_ptr),
            ),
        ]:
            exception = {}
            error = make_interface(exception)(*args, None)
            if error != 1 or not isinstance(exception.get("exc"), ValueError):
                print(
                    f" test_callback_interfaces failed: Writing into the input of the "
                    f"{name} not reported as error."
                )
                test_passed = False
            if np.any(vector != 0.0):
                print(f" test_callback_interfaces failed: Input of the {name} changed.")
                test_passed = False
                vector[:] = 0.0

        # the orbital update interface leaves the Hessian linear transformation unset
        # without an error when the orbital update does not provide one, so that the
        # solver reports the missing one, and reports one that is not callable as an
        # error of the orbital update
        for returned, expected_error, case in [
            (None, 0, "without Hessian linear transformation"),
            (
                np.zeros(n_param),
                1,
                "with Hessian linear transformation that is not callable",
            ),
        ]:
            exception = {}
            hess_x_funptr = pointer(hess_x_interface_type())
            error = UpdateOrbsInterface(
                lambda kappa, grad, h_diag: (0.0, returned), n_param, exception
            ).call(vector_ptr, pointer(c_real()), vector_ptr, vector_ptr, hess_x_funptr)
            if error != expected_error:
                print(
                    f" test_callback_interfaces failed: Error reported wrong for "
                    f"orbital update {case}."
                )
                test_passed = False
            if expected_error and not isinstance(exception.get("exc"), TypeError):
                print(
                    f" test_callback_interfaces failed: Exception not kept for "
                    f"orbital update {case}."
                )
                test_passed = False
            if hess_x_funptr[0]:
                print(
                    f" test_callback_interfaces failed: Hessian linear transformation "
                    f"set for orbital update {case}."
                )
                test_passed = False

        # the Hessian linear transformation an orbital update returns reports an
        # exception it raises as an error and keeps it with those of the orbital update
        exception = {}
        hess_x_funptr = pointer(hess_x_interface_type())
        UpdateOrbsInterface(
            lambda kappa, grad, h_diag: (0.0, raising(ValueError)), n_param, exception
        ).call(vector_ptr, pointer(c_real()), vector_ptr, vector_ptr, hess_x_funptr)
        if hess_x_funptr[0](vector_ptr, output_ptr, None) != 1 or not isinstance(
            exception.get("exc"), ValueError
        ):
            print(
                " test_callback_interfaces failed: Exception raised by Hessian linear "
                "transformation returned by orbital update not reported."
            )
            test_passed = False

        # the convergence check interface passes on that the optimization has not
        # converged
        converged = c_bool(True)
        if ConvCheckInterface(lambda: False, {}).call(pointer(converged)) != 0:
            print(
                " test_callback_interfaces failed: Error reported for convergence "
                "check."
            )
            test_passed = False
        if converged.value:
            print(
                " test_callback_interfaces failed: Result of convergence check not "
                "passed on."
            )
            test_passed = False

        # the logging interface keeps the first of several exceptions
        failures = iter([ValueError("first"), ValueError("second")])

        def failing_logger(message):
            raise next(failures)

        exception = {}
        logger_interface = LoggerInterface(failing_logger, exception)
        logger_interface(b"test", None)
        logger_interface(b"test", None)
        if str(exception.get("logger")) != "first":
            print(
                " test_callback_interfaces failed: First exception raised by "
                "logging function not kept."
            )
            test_passed = False

        # an exception raised by an interface itself, as for a pointer it can not
        # write the result through, is reported as an error and kept
        exception = {}
        if ObjFuncInterface(lambda kappa: 0.0, n_param, exception)(
            vector_ptr, None, None
        ) != 1 or not isinstance(exception.get("exc"), TypeError):
            print(
                " test_callback_interfaces failed: Exception raised by interface not "
                "reported as error."
            )
            test_passed = False

        # an interface marks a callback function as innermost while it runs and
        # unmarks it afterwards, also when it raises
        try:
            for error_type, case in [(None, "returns"), (ValueError, "raises")]:
                interruption_state.stack = (False,)
                marks = []

                def marking(kappa):
                    marks.append(interruption_state.stack)
                    if error_type is not None:
                        raise error_type("objective function failure")
                    return 0.0

                ObjFuncInterface(marking, n_param, {})(
                    vector_ptr, pointer(c_real()), None
                )
                if marks != [(False, True)] or interruption_state.stack != (False,):
                    print(
                        " test_callback_interfaces failed: Callback function not "
                        f"marked as innermost while it runs when it {case}."
                    )
                    test_passed = False
        finally:
            interruption_state.stack = ()

        # the logging interface records an interruption, which the next callback
        # function with an error code reports without being called
        try:
            interruption = KeyboardInterrupt()

            def interrupted_logger(message):
                raise interruption

            exception = {}
            LoggerInterface(interrupted_logger, exception)(b"test", None)
            LoggerInterface(lambda message: None, exception)(b"test", None)
            if interruption_state.pending is not interruption or "logger" in exception:
                print(
                    " test_callback_interfaces failed: Interruption of logging "
                    "function not recorded."
                )
                test_passed = False
            obj_func_calls = []
            error = ObjFuncInterface(obj_func_calls.append, n_param, exception)(
                vector_ptr, pointer(c_real()), None
            )
            if error != 1 or exception.get("exc") is not interruption or obj_func_calls:
                print(
                    " test_callback_interfaces failed: Recorded interruption not "
                    "reported as error by next callback function."
                )
                test_passed = False
            if interruption_state.pending is not None:
                print(
                    " test_callback_interfaces failed: Recorded interruption kept "
                    "after it was reported."
                )
                test_passed = False
        finally:
            interruption_state.pending = None

        self.assertTrue(test_passed, "test_callback_interfaces failed")
        print(" test_callback_interfaces PASSED")

    def test_raise_on_failure(self):
        """
        this function tests that an error returned by the library is raised with its
        code and the exception of a failed callback function as its cause, that a
        failed logging function is reported as a warning with its exception as the
        cause both with and without an error, that nothing is reported without either
        and that an interruption of a callback function is raised as it is
        """
        test_passed = True
        callback_failure = ValueError("callback function failure")
        logger_failure = ValueError("logging failure")

        for error, exception, case in [
            (7, {}, "for an error alone"),
            (
                7,
                {"exc": callback_failure},
                "for an error of a failed callback function",
            ),
            (7, {"logger": logger_failure}, "for an error and a failed logger"),
            (0, {"logger": logger_failure}, "for a failed logger alone"),
            (0, {}, "without failure"),
        ]:
            raised = None
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                try:
                    raise_on_failure(error, exception, "solver")
                except RuntimeError as e:
                    raised = e

            # an error is raised with its code and the exception of a failed callback
            # function as its cause, never the one of a failed logging function
            if error and (raised is None or "(code 7)" not in str(raised)):
                print(
                    f" test_raise_on_failure failed: Error not raised with its code "
                    f"{case}."
                )
                test_passed = False
            if raised is not None and raised.__cause__ is not exception.get("exc"):
                print(
                    f" test_raise_on_failure failed: Error raised with wrong cause "
                    f"{case}."
                )
                test_passed = False
            if not error and raised is not None:
                print(f" test_raise_on_failure failed: Raised without error {case}.")
                test_passed = False

            # a failed logging function is reported as a warning with its exception as
            # the cause and nothing else is warned about
            warned_causes = [
                w.message.__cause__
                for w in caught
                if issubclass(w.category, RuntimeWarning)
            ]
            if "logger" in exception and warned_causes != [logger_failure]:
                print(
                    f" test_raise_on_failure failed: Failed logging function not "
                    f"reported as a warning {case}."
                )
                test_passed = False
            if "logger" not in exception and caught:
                print(
                    f" test_raise_on_failure failed: Warned without failed logger {case}."
                )
                test_passed = False

        # an interruption of a callback function is raised as it is, rather than as an
        # error of the library
        interruption = KeyboardInterrupt()
        raised = None
        try:
            raise_on_failure(7, {"exc": interruption}, "solver")
        except BaseException as e:
            raised = e
        if raised is not interruption:
            print(
                " test_raise_on_failure failed: Interruption of callback function not "
                "raised as it is."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_raise_on_failure failed")
        print(" test_raise_on_failure PASSED")

    def test_call_library(self):
        """
        this function tests that a call into the library passes on the arguments and
        the error code and marks the library as innermost, that it replaces Python's
        default handler of SIGINT in the main thread only and restores it afterwards,
        also when the call raises, and that it raises an interruption that arrived
        during the call
        """
        test_passed = True
        seen = []

        def library(*args):
            seen.append(
                (args, interruption_state.stack, signal.getsignal(signal.SIGINT))
            )
            return 7

        def check_restored(case):
            if signal.getsignal(signal.SIGINT) is not signal.default_int_handler:
                print(f" test_call_library failed: Handler not restored {case}.")
                return False
            if interruption_state.stack != ():
                print(f" test_call_library failed: Marks not restored {case}.")
                return False
            return True

        original_handler = signal.signal(signal.SIGINT, signal.default_int_handler)
        try:
            # the outermost call from the main thread defers interruptions
            seen.clear()
            if call_library(library, 1, "a") != 7 or seen != [
                ((1, "a"), (False,), defer_interruption)
            ]:
                print(
                    " test_call_library failed: Outermost call from the main thread "
                    "not made correctly."
                )
                test_passed = False
            test_passed &= check_restored("after the outermost call")

            # a nested call keeps the handler and nests its mark
            seen.clear()
            call_library(lambda: call_library(library))
            if seen != [((), (False, False), defer_interruption)]:
                print(" test_call_library failed: Nested call not made correctly.")
                test_passed = False
            test_passed &= check_restored("after a nested call")

            # a call from another thread neither replaces the handler nor marks the
            # main thread
            seen.clear()
            thread = threading.Thread(target=call_library, args=(library,))
            thread.start()
            thread.join()
            if seen != [((), (False,), signal.default_int_handler)]:
                print(
                    " test_call_library failed: Call from another thread not made "
                    "correctly."
                )
                test_passed = False
            test_passed &= check_restored("after a call from another thread")

            # a call raises the exception of the library and restores the handler
            def failing_library():
                raise ValueError("library failure")

            try:
                call_library(failing_library)
                print(" test_call_library failed: Exception of library not raised.")
                test_passed = False
            except ValueError:
                pass
            test_passed &= check_restored("after a call that raises")

            # an interruption the library received is raised when the call returns
            interruption = KeyboardInterrupt()

            def interrupted_library():
                interruption_state.pending = interruption
                return 0

            raised = None
            try:
                call_library(interrupted_library)
            except KeyboardInterrupt as e:
                raised = e
            if raised is not interruption or interruption_state.pending is not None:
                print(
                    " test_call_library failed: Interruption received by the library "
                    "not raised when the call returns."
                )
                test_passed = False
            test_passed &= check_restored("after an interruption")

            # an interruption that arrives while the handler is restored, which Python
            # handles while it switches the handler, is raised when the call returns,
            # after the handler is restored
            real_signal = signal.signal

            def signal_handling_pending(signum, handler):
                if handler is signal.default_int_handler:
                    signal.getsignal(signum)(signum, None)
                return real_signal(signum, handler)

            raised = None
            with patch.object(signal, "signal", signal_handling_pending):
                try:
                    call_library(library)
                except KeyboardInterrupt as e:
                    raised = e
            if raised is None or interruption_state.pending is not None:
                print(
                    " test_call_library failed: Interruption that arrives while the "
                    "handler is restored not raised when the call returns."
                )
                test_passed = False
            test_passed &= check_restored("after an interruption while it is restored")

            # another handler is kept
            def custom_handler(signum, frame):
                pass

            signal.signal(signal.SIGINT, custom_handler)
            seen.clear()
            call_library(library)
            if seen[0][2] is not custom_handler:
                print(" test_call_library failed: Other handler replaced.")
                test_passed = False
            if signal.getsignal(signal.SIGINT) is not custom_handler:
                print(" test_call_library failed: Other handler not kept.")
                test_passed = False
        finally:
            signal.signal(signal.SIGINT, original_handler)
            interruption_state.stack = ()
            interruption_state.pending = None

        self.assertTrue(test_passed, "test_call_library failed")
        print(" test_call_library PASSED")

    def test_check_struct_sizes(self):
        """
        this function tests that the check of the structure sizes passes for a library
        whose structures have the sizes of the Python structures and raises an import
        error for one whose structure differs in size or which exports no sizes
        """
        test_passed = True
        sizes = {
            "solver_settings": sizeof(SolverSettings.c_struct),
            "stability_settings": sizeof(StabilitySettings.c_struct),
            "solver_state": sizeof(SolverState.c_struct),
            "stability_state": sizeof(StabilityState.c_struct),
        }

        def mock_library(sizes):
            """
            this function returns a mock library which exports the given sizes
            """
            return SimpleNamespace(
                _name="mock library",
                **{
                    f"{name}_size": (lambda size=size: size)
                    for name, size in sizes.items()
                },
            )

        # the check passes for the sizes of the Python structures
        try:
            check_struct_sizes(mock_library(sizes))
        except ImportError:
            print(" test_check_struct_sizes failed: Raised for matching sizes.")
            test_passed = False

        # the check raises for a structure of another size and for a missing size
        for name in sizes:
            for wrong_sizes, case in [
                ({**sizes, name: sizes[name] + 8}, f"for a differing {name} size"),
                (
                    {key: size for key, size in sizes.items() if key != name},
                    f"for a missing {name} size",
                ),
            ]:
                try:
                    check_struct_sizes(mock_library(wrong_sizes))
                    print(f" test_check_struct_sizes failed: Not raised {case}.")
                    test_passed = False
                except ImportError:
                    pass

        self.assertTrue(test_passed, "test_check_struct_sizes failed")
        print(" test_check_struct_sizes PASSED")

    @patch.object(SolverSettings, "init_c_struct", lib.mock_init_solver_settings)
    def test_solver_settings(self):
        """
        this function ensures the SolverSettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = SolverSettings()
        test_passed = self._settings_hold_reference(
            "test_solver_settings", settings, ""
        )
        test_passed &= self._settings_safeguarded("test_solver_settings", settings)
        test_passed &= self._callbacks_wired("test_solver_settings", SolverSettings)

        dummy_error_code = 42

        def dummy_precond():
            return dummy_error_code

        # an optional callback function is set in the C struct of a call and its
        # interface is kept in the interfaces of the call
        settings_c = settings.c_struct()
        interfaces = []
        settings.set_optional_callback(
            settings_c,
            interfaces,
            "precond",
            dummy_precond,
            lambda x: x,
            CFUNCTYPE(c_int),
        )

        c_ptr = settings_c.precond
        c_interface = interfaces[-1] if interfaces else None

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

        # an unset optional callback function leaves its C function pointer null
        settings.set_optional_callback(
            settings_c, interfaces, "precond", None, lambda x: x, CFUNCTYPE(c_int)
        )
        if settings_c.precond is not None or interfaces[-1] is not None:
            print(
                " test_solver_settings failed: Unset optional callback does not "
                "leave a null function pointer."
            )
            test_passed = False

        # a keyword given as bytes, as read from a C struct, is stored as well
        settings.subsystem_solver = b"davidson_ah"
        if settings.subsystem_solver != "davidson_ah":
            print(
                " test_solver_settings failed: Keyword given as bytes not stored "
                "correctly."
            )
            test_passed = False

        self.assertTrue(test_passed, "test_solver_settings failed")
        print(" test_solver_settings PASSED")

    @patch.object(StabilitySettings, "init_c_struct", lib.mock_init_stability_settings)
    def test_stability_settings(self):
        """
        this function ensures the StabilitySettings object is properly initialized and
        synchronized with the underlying C struct
        """
        settings = StabilitySettings()
        test_passed = self._settings_hold_reference(
            "test_stability_settings", settings, "stability_settings."
        )
        test_passed &= self._settings_safeguarded("test_stability_settings", settings)
        test_passed &= self._callbacks_wired(
            "test_stability_settings", StabilitySettings
        )
        self.assertTrue(test_passed, "test_stability_settings failed")
        print(" test_stability_settings PASSED")

    def test_solver_state(self):
        """
        this function ensures the SolverState object starts from the default values
        and reads every field of the underlying C struct, including those of the
        nested stability check state, which it does not let be set
        """
        state = SolverState()
        test_passed = self._state_holds_defaults("test_solver_state", state)
        lib.get_reference_solver_state_values(byref(state.state_c))
        test_passed &= self._state_holds_reference("test_solver_state", state, "")
        test_passed &= self._state_copies_kept("test_solver_state", state, "")
        try:
            state.n_update_orbs = 1
            print(" test_solver_state failed: Field could be set.")
            test_passed = False
        except AttributeError:
            pass
        self.assertTrue(test_passed, "test_solver_state failed")
        print(" test_solver_state PASSED")

    def test_stability_state(self):
        """
        this function ensures the StabilityState object starts from the default values
        and reads every field of the underlying C struct, which it does not let be set
        """
        state = StabilityState()
        test_passed = self._state_holds_defaults("test_stability_state", state)
        lib.get_reference_stability_state_values(byref(state.state_c))
        test_passed &= self._state_holds_reference(
            "test_stability_state", state, "stability_state."
        )
        test_passed &= self._state_copies_kept(
            "test_stability_state", state, "stability_state."
        )
        try:
            state.n_hess_x = 1
            print(" test_stability_state failed: Field could be set.")
            test_passed = False
        except AttributeError:
            pass
        self.assertTrue(test_passed, "test_stability_state failed")
        print(" test_stability_state PASSED")


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

    def test_state_layout(self):
        """
        this function checks that every field of the Python state structures,
        including the nested stability check state, is read back under its own name
        from a state whose fields Fortran sets one by one to the reference values
        """
        test_passed = True

        # check that every field is read back under its own name
        state_c = SolverState.c_struct()
        lib.get_reference_solver_state_values(byref(state_c))
        for struct, prefix in [
            (state_c, ""),
            (state_c.stability_state, "stability_state."),
        ]:
            for name, field_type in struct._fields_:
                if field_type == StabilityState.c_struct:
                    continue
                if getattr(struct, name) != reference_state(prefix + name, field_type):
                    print(
                        f" test_state_layout failed: Field {prefix + name} misplaced."
                    )
                    test_passed = False
        self.assertTrue(test_passed, "test_state_layout failed")
        print(" test_state_layout PASSED")

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
            "n_update_orbs_calls": 0,
            "state_not_current": False,
            "nested_n_hess_x_seen": None,
            "residual_rms_seen": None,
        }

        # the solver state record, which the callbacks read during the solve
        solver_state = SolverState()

        def update_orbs(delta_vars, grad, h_diag):
            # the state is copied before every callback, so it counts the updates
            # before this one
            if solver_state.n_update_orbs != state["n_update_orbs_calls"]:
                state["state_not_current"] = True
            state["n_update_orbs_calls"] += 1
            state["curr"] = state["curr"] + delta_vars
            x = state["curr"]
            state["hess"] = self._hess(x)
            grad[:] = self._grad(x)
            h_diag[:] = np.diag(state["hess"])

            def hess_x(v, hv):
                # the internal stability check calls this function too, the last time
                # after all its other Hessian linear transformations
                state["nested_n_hess_x_seen"] = solver_state.stability_state.n_hess_x
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
            # report convergence by the gradient criterion with the default threshold
            # it replaces, so the solve is unaffected
            state["called"].add("conv_check")
            return (
                np.sqrt(np.mean(self._grad(state["curr"]) ** 2))
                < SolverSettings().conv_tol
            )

        # extra trial vector functions that contribute no direction, so that the solve
        # is unaffected
        def recording_get_extra_trial_vectors(name):
            def get_extra_trial_vectors(trial_vectors):
                state["called"].add(name)
                trial_vectors[:] = 0.0

            return get_extra_trial_vectors

        # convergence check of the internal stability check, which reads the RMS
        # residual of the lowest eigenpair from the nested state, records it and
        # accepts a converged residual
        def stability_conv_check():
            state["called"].add("stability_conv_check")
            state["residual_rms_seen"] = solver_state.stability_state.residual_rms
            return solver_state.stability_state.residual_rms < 1e-9

        settings = SolverSettings()
        settings.precond = recording_precond("precond")
        settings.project = recording_callback("project")
        settings.get_extra_trial_vectors = recording_get_extra_trial_vectors(
            "get_extra_trial_vectors"
        )
        settings.conv_check = conv_check
        settings.logger = recording_callback("logger")
        settings.stability = True
        settings.stability_settings.precond = recording_precond("stability_precond")
        settings.stability_settings.project = recording_callback("stability_project")
        settings.stability_settings.get_extra_trial_vectors = (
            recording_get_extra_trial_vectors("stability_get_extra_trial_vectors")
        )
        settings.stability_settings.conv_check = stability_conv_check
        settings.stability_settings.logger = recording_callback("stability_logger")
        settings.verbose = 3  # ensure the loggers are exercised

        # the maximum precision flag starts opposite to the one the converging solve
        # returns, so that its write-back is detected
        solver_state.state_c.max_precision_reached = True

        solver(update_orbs, obj_func, self.n_param, settings, state=solver_state)
        if solver_state.max_precision_reached:
            print(
                " test_solver_py failed: Maximum precision reached flag was not "
                "returned."
            )
            test_passed = False
        for name, description in [
            ("precond", "Preconditioner"),
            ("project", "Projection"),
            ("get_extra_trial_vectors", "Extra trial vector function"),
            ("conv_check", "Convergence check"),
            ("logger", "Logger"),
        ]:
            if name not in state["called"]:
                print(f" test_solver_py failed: {description} was not called.")
                test_passed = False
        if solver_state.n_update_orbs <= 0:
            print(" test_solver_py failed: Orbital update counter was not populated.")
            test_passed = False
        if solver_state.n_obj_func <= 0:
            print(
                " test_solver_py failed: Objective function counter was not populated."
            )
            test_passed = False
        if solver_state.n_hess_x <= 0:
            print(
                " test_solver_py failed: Hessian linear transformation counter was not "
                "populated."
            )
            test_passed = False
        if solver_state.stability_state.n_hess_x <= 0:
            print(
                " test_solver_py failed: Hessian linear transformation counter of the "
                "internal stability check was not populated."
            )
            test_passed = False

        # every orbital update read the updates before it, the last Hessian linear
        # transformation of the internal stability check the ones before it and its
        # last convergence check the residual of its last iteration
        if (
            state["state_not_current"]
            or state["n_update_orbs_calls"] != solver_state.n_update_orbs
            or state["nested_n_hess_x_seen"]
            != solver_state.stability_state.n_hess_x - 1
            or state["residual_rms_seen"] != solver_state.stability_state.residual_rms
        ):
            print(" test_solver_py failed: State was not current in the callbacks.")
            test_passed = False
        for name, description in [
            ("stability_precond", "Preconditioner"),
            ("stability_project", "Projection"),
            ("stability_get_extra_trial_vectors", "Extra trial vector function"),
            ("stability_conv_check", "Convergence check"),
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

        # an extra trial vector function that contributes no direction and a
        # convergence check that reads the RMS residual of the lowest eigenpair from the
        # state, records it and accepts a converged residual
        def get_extra_trial_vectors(trial_vectors):
            called.add("get_extra_trial_vectors")
            trial_vectors[:] = 0.0

        residual_rms_seen = []

        def conv_check():
            called.add("conv_check")
            residual_rms_seen.append(stability_state.residual_rms)
            return stability_state.residual_rms < 1e-9

        settings = StabilitySettings()
        settings.precond = precond
        settings.project = project
        settings.get_extra_trial_vectors = get_extra_trial_vectors
        settings.conv_check = conv_check
        settings.logger = logger
        settings.verbose = 3  # ensure logger is exercised

        # at a minimum, expect stable
        H = self._hess(self.minimum1)
        h_diag = np.diag(H).copy()

        def hess_x(v, hv):
            hv[:] = H @ v

        # the stability check state record, which the Hessian linear transformation
        # reads during the first call
        stability_state = StabilityState()
        n_hess_x_seen = []

        def watching_hess_x(v, hv):
            n_hess_x_seen.append(stability_state.n_hess_x)
            hess_x(v, hv)

        kappa = np.zeros(self.n_param)
        stable = stability_check(
            h_diag,
            watching_hess_x,
            self.n_param,
            settings,
            kappa=kappa,
            state=stability_state,
        )
        if not stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of minimum."
            )
            test_passed = False
        if stability_state.iteration <= 0:
            print(" test_stability_check_py failed: State was not returned.")
            test_passed = False

        # the last Hessian linear transformation read the ones before it and the last
        # convergence check the residual of the last iteration
        if (
            not n_hess_x_seen
            or n_hess_x_seen[-1] != stability_state.n_hess_x - 1
            or not residual_rms_seen
            or residual_rms_seen[-1] != stability_state.residual_rms
        ):
            print(
                " test_stability_check_py failed: State was not current in the "
                "callbacks."
            )
            test_passed = False
        for name, description in [
            ("precond", "Preconditioner"),
            ("project", "Projection"),
            ("get_extra_trial_vectors", "Extra trial vector function"),
            ("conv_check", "Convergence check"),
            ("logger", "Logger"),
        ]:
            if name not in called:
                print(f" test_stability_check_py failed: {description} was not called.")
                test_passed = False
        if stability_state.n_hess_x <= 0:
            print(
                " test_stability_check_py failed: Hessian linear transformation "
                "counter was not populated."
            )
            test_passed = False
        if stability_state.min_eigval <= 0.0:
            print(
                " test_stability_check_py failed: State does not hold a positive "
                "lowest eigenvalue for minimum."
            )
            test_passed = False

        # also exercise the call without a direction at the minimum, where the
        # stability flag differs from the False the wrapper starts from, with the
        # state, which the convergence check reads
        stable = stability_check(
            h_diag, hess_x, self.n_param, settings, state=stability_state
        )
        if not stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of minimum when not passing direction."
            )
            test_passed = False

        # at a saddle, the stability flag and the sign of the lowest eigenvalue have
        # to differ from the ones returned at the minimum and the returned direction
        # has to be a normalized direction of negative curvature
        H = self._hess(self.saddle_point)
        h_diag = np.diag(H).copy()
        stable = stability_check(
            h_diag, hess_x, self.n_param, settings, kappa=kappa, state=stability_state
        )
        if stable:
            print(
                " test_stability_check_py failed: Stability incorrectly classifies "
                "stability of saddle point."
            )
            test_passed = False
        if stability_state.min_eigval >= 0.0:
            print(
                " test_stability_check_py failed: State does not hold a negative "
                "lowest eigenvalue for saddle point."
            )
            test_passed = False
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
