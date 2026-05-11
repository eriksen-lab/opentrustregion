# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

OpenTrustRegion is a Fortran library implementing second-order trust region orbital optimization for quantum chemistry, exposed via Fortran, C, and Python (ctypes) interfaces. The same compiled shared/static library is consumed by all three.

## Project priorities

**Clarity over performance.** The library is not on the performance hot path of a quantum-chemistry calculation — the host program's integral transforms and Hessian-vector products dominate. Prefer code that is short and obviously correct over code that is fast. Don't propose performance refactors (amortized buffer growth, pooling, micro-optimizations) unless there is evidence the affected code is hot for a real workload.

**Fortran/C/Python interfaces must always be consistent.** The same callback signatures, settings fields, and default values are described in six places that must agree:

1. Fortran abstract interfaces and `solver_settings_type` / `stability_settings_type` in `src/opentrustregion.f90`
2. C abstract interfaces and `bind(C)` `solver_settings_type_c` / `stability_settings_type_c in `src/c_interface.f90` 
3. C struct layouts and typedefs in `include/opentrustregion.h`
4. `SolverSettingsC` / `StabilitySettingsC` ctypes `_fields_` and `CFUNCTYPE` declarations in `pyopentrustregion/python_interface.py`
5. `default_solver_settings` / `default_stability_settings` (Fortran) ↔ `solver_settings_init` / `stability_settings_init` (C) ↔ the Python `Settings` wrapper defaults
6. The argument lists and snippets in `README.md`

Any change to a callback signature, a settings field (add / remove / reorder), or a default value must land in all six locations in the same PR. Additionally, the error-origin codes (`error_obj_func` etc. in `opentrustregion.f90`) must be synchronized with the README error-code table.

## Build & test

```sh
# Fortran/C build
mkdir build && cd build
cmake ..              # add -DBUILD_SHARED_LIBS=ON for shared
cmake --build .

# Python install (invokes CMake under the hood via setup.py)
pip install .
pip install -e .       # editable
```

Run the full test suite (Python ctypes driver that calls into `libotrtestsuite`, which contains the Fortran unit + system tests):

```sh
python3 -m pyopentrustregion.testsuite                  # from an installed/editable build
python3 pyopentrustregion/testsuite.py                  # from the source tree against ./build
```

Run a single test class or method using stdlib `unittest`:

```sh
python3 -m unittest pyopentrustregion.testsuite.SystemTests
python3 -m unittest pyopentrustregion.testsuite.OpenTrustRegionUnitTests.test_solver
```

The Python suite runs Fortran- and C-side tests (via symbols dynamically loaded from `libotrtestsuite`) and pure-Python wrapper tests. System tests need `pyopentrustregion/test_data/*.bin`.

### CMake options that matter

- `INTEGER_SIZE` (`4` or `8`): library integer width. If unset, CMake autodetects by trying 32-bit BLAS/LAPACK first, falling back to 64-bit. Output library is named `libopentrustregion_32.*` or `libopentrustregion_64.*`. Defining `USE_ILP64` (set automatically when `INTEGER_SIZE=8`) switches both Fortran `ip` and the C `c_ip` typedefs to 64-bit and remaps BLAS/LAPACK symbol names to their `_64` variants when `check_fortran_function_exists` finds them.
- `BLAS_LIBRARIES` / `LAPACK_LIBRARIES`: required to be set together with `INTEGER_SIZE` if overriding autodetection — passing one without the other is a fatal error.
- `OpenTrustRegion_BUILD_TESTING` (default ON when top-level): builds `libotrtestsuite`, the shared library Python loads to drive the Fortran tests.
- `CONDA_BUILD=1` env var: the Python `setup.py` skips the embedded CMake invocation (the conda recipe builds the C library separately).

## Architecture

### Source layout

- `src/opentrustregion.f90` — the entire numerical core: `solver`, `stability_check`, the Davidson / Jacobi-Davidson / truncated-CG subsystem solvers, settings derived types, callback abstract interfaces, error-code constants. ~2.6k lines, single module.
- `src/c_interface.f90` — `bind(C)` wrapper module. Stores Fortran procedure pointers to C callbacks at module scope (`update_orbs_before_wrapping`, etc.) and adapts C-style `(*)` arrays + return-code functions into the Fortran-style `(:)` arrays + `intent(out) :: error` subroutines that the core expects.
- `include/opentrustregion.h` — C header mirroring the Fortran settings types (`solver_settings_type`, `stability_settings_type`) as C structs plus `*_init()` helpers and `solver` / `stability_check` prototypes.
- `pyopentrustregion/python_interface.py` — ctypes wrapper. Defines `SolverSettings` / `StabilitySettings` as `ctypes.Structure` mirrors of the C structs, wraps Python callbacks with `CFUNCTYPE`, and converts the integer error codes returned by the C entry points into `RuntimeError`.
- `tests/` — three layers:
  - `opentrustregion_unit_tests.f90` / `c_interface_unit_tests.f90` — unit tests against the Fortran and C interfaces respectively, both built into `libotrtestsuite`.
  - `opentrustregion_system_tests.f90` — system tests using reference binary inputs from `pyopentrustregion/test_data/`, built into `libotrtestsuite`.
  - `c_system_tests.c` — system tests through the C interface.
  - `*_mock.f90` — mock callbacks used by both unit test suites.
  - `test_reference.f90` — shared tolerance constants and reference values.

### Three-language interface chain

Python → C → Fortran is one call path:
`solver(...)` in `python_interface.py` builds a `SolverSettings` ctypes struct, wraps the user's Python callbacks with `CFUNCTYPE`, and calls into `libopentrustregion`'s C entry point. That entry point lives in `src/c_interface.f90`, where the Fortran wrapper stores the C function pointers in module-level `procedure(..)_before_wrapping` slots, then invokes `standard_solver` (the real `solver` from the `opentrustregion` module) with Fortran-shaped wrapper procedures that internally call back through those module-level pointers. **Consequence:** the C interface is not re-entrant across threads — module-level pointer state is shared. When changing any signature, all three layers (`opentrustregion.f90` abstract interface → `c_interface.f90` C wrapper + abstract interface → `opentrustregion.h` typedef → `python_interface.py` `CFUNCTYPE` + `Structure`) must be kept in lockstep.

### Settings types

Both `solver_settings_type` and `stability_settings_type` extend an abstract `settings_type` and have an `init` type-bound procedure plus an `initialized` flag. The solver checks `settings%initialized` on entry and calls `init` if needed; the C and Python wrappers pre-populate defaults via their respective `*_init()` / `__init__` methods so that flag is `.true.` by the time they reach Fortran. Default values are defined once in `default_solver_settings` / `default_stability_settings` parameters in `opentrustregion.f90` and replicated in `solver_settings_init()` (C) and the Python `Structure` defaults — keep them in sync.

### Error codes

Errors are encoded as `OOEE` integers (origin × 100 + specific code, see README). The Fortran core uses `add_error_origin` to tag a non-zero error coming back from a callback with the appropriate origin (e.g. `error_update_orbs = 1200`). The C entry points return this integer directly; the Python wrappers check it and raise `RuntimeError` with a decoded origin/code. New origins must be added as a parameter in `opentrustregion.f90` and decoded in `python_interface.py`.

### Preprocessing

All Fortran sources are compiled with `Fortran_PREPROCESS ON`. Integer kind selection (`int32` vs `int64`) and the BLAS/LAPACK 64-bit symbol remapping (`ddot=ddot_64`, etc.) happen via `#ifdef USE_ILP64` and `add_compile_definitions` in CMake — do not hardcode integer kinds.
