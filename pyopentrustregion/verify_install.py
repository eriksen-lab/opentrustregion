# Copyright (C) 2025- Jonas Greiner
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

"""
Verify that an installed OpenTrustRegion works, via the public Python API.

Run as "python -m pyopentrustregion.verify_install".

The wheels published on PyPI carry a vendored BLAS/LAPACK and Fortran runtime but omit
the Fortran testsuite library, so pyopentrustregion.testsuite cannot run against them.
This drives the solver through the public API instead and checks that it returns
correct numbers, which is what catches an installation whose BLAS is mislinked or built
for the wrong integer width. Source and conda installs can run the full testsuite as
well.
"""

from __future__ import annotations

import numpy as np
from typing import TYPE_CHECKING

from pyopentrustregion import SolverSettings, solver

if TYPE_CHECKING:
    from typing import Tuple, Callable, Dict

# minimize sum(x**2) starting away from its minimum at x = 0
N_PARAM: int = 6
state: Dict[str, np.ndarray] = {"curr": np.array([0.20, 0.15, 0.48, 0.28, 0.31, 0.66])}


def obj_func(delta_vars: np.ndarray) -> float:
    x = state["curr"] + delta_vars
    return float(np.dot(x, x))


def update_orbs(
    delta_vars: np.ndarray, grad: np.ndarray, h_diag: np.ndarray
) -> Tuple[float, Callable[[np.ndarray, np.ndarray], None]]:
    state["curr"] = state["curr"] + delta_vars
    x = state["curr"]
    grad[:] = 2.0 * x
    h_diag[:] = 2.0

    def hess_x(vector: np.ndarray, hess_vector: np.ndarray) -> None:
        hess_vector[:] = 2.0 * vector

    return float(np.dot(x, x)), hess_x


def main() -> None:
    # run the solver
    settings = SolverSettings()
    settings.conv_tol = 1e-8
    settings.n_macro = 100
    settings.verbose = 0
    solver(obj_func, update_orbs, N_PARAM, settings)

    # verify the solver actually reached the known minimum
    x = state["curr"]
    if not np.allclose(x, 0.0, atol=1e-6):
        raise SystemExit(
            f"OpenTrustRegion installation check FAILED: solver returned x = {x}"
        )

    print("OpenTrustRegion installation check passed.")


if __name__ == "__main__":
    main()
