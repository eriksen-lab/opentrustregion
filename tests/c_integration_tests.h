// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* Declarations of the C integration tests: the capture of the standard error, which
 * c_integration_tests.c provides, and the warnings the library prints for settings
 * that were not initialized. */

#ifndef C_INTEGRATION_TESTS_H
#define C_INTEGRATION_TESTS_H

#include <stdbool.h>

/* the warnings the library prints for settings and nested stability check settings
 * that were not initialized, provided by test_reference.f90 */
extern const char test_settings_uninitialized_warning[];
extern const char test_stability_settings_uninitialized_warning[];

/* run the calls between these two functions with the standard error redirected to a
 * temporary file, flushing what C and the Fortran runtime buffered before and after
 * them, the second returns whether the given warning was written to it, or with NULL
 * whether nothing was written to it */
void start_stderr_capture(void);
bool stop_stderr_capture(const char *warning);

#endif
