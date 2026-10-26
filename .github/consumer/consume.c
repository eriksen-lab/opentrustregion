// Copyright (C) 2025- Jonas Greiner
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/* ------------------------------------------------------------------
 * A C program built against an installed OpenTrustRegion, which runs a stability
 * check on a diagonal Hessian with positive eigenvalues and fails unless it finds the
 * point stable without an error
 * ------------------------------------------------------------------ */

#include "opentrustregion.h"
#include <stdio.h>

/* Hessian linear transformation of the diagonal Hessian diag(1, 2, 3) */
static c_int hess_x(const c_real *x, c_real *hx, void *context) {
  (void)context;
  for (int i = 0; i < 3; i++)
    hx[i] = (i + 1) * x[i];
  return 0;
}

int main(void) {
  c_real h_diag[3] = {1.0, 2.0, 3.0};
  c_bool stable;
  c_int error = stability_check(h_diag, hess_x, 3, &stable, NULL, NULL, NULL);
  printf("C program: error %d, stable %d\n", (int)error, (int)stable);
  return error != 0 || !stable;
}
