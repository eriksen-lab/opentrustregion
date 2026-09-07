# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

OpenTrustRegion is a Fortran library implementing a second-order trust-region optimizer, exposed via Fortran, C, and Python (ctypes) interfaces — the same compiled shared/static library is consumed by all three.

## General conventions

**Formatting.** Fortran lines max 88 columns, continued with a trailing `&` aligned to the opening parenthesis of the call. Every procedure opens with a `!`-delimited comment block describing what it does; each logical step inside gets a lowercase `!` comment. Unit tests: assume success (`test_<name> = .true.`), then one `if (...) then` / `write (stderr, *) "test_<name> failed: ..."` / `test_<name> = .false.` block per assertion, no early returns except where a later assertion would crash.

**Python formatting.** All Python (`pyopentrustregion/`, `setup.py`) is `black`-formatted, default settings. Run `black pyopentrustregion setup.py` before considering Python changes done.

**Clarity over performance.** The library is not on the hot path of a quantum-chemistry calculation — the host program's integral transforms and Hessian-vector products dominate. Prefer short, obviously-correct code over fast code. Don't propose performance refactors (buffer growth, pooling, micro-optimizations) without evidence the affected code is hot for a real workload.

**Norms: BLAS in production, intrinsics in tests.** Production Fortran (`src/`) uses BLAS (`dnrm2`, `ddot`) for norms/dot products. Test code (`tests/`) uses intrinsics (`norm2`, `dot_product`, `matmul`, `sum`, `transpose`) instead. BLAS/LAPACK is allowed in tests only where no intrinsic exists (`dsyev`, `dgeev`, `zheev`) — don't hand-roll linear algebra to avoid the dependency.

### Unit tests must not depend on untested code

A unit test may only call the routine it is testing — never another production routine to build its inputs or expected values, or the test silently inherits that routine's correctness.

- Construct inputs directly in the test, or with LAPACK, not by calling the production routine that would normally produce them.
- When a routine under test internally calls another module's routines, reimplement those as `ref_*` helpers so the expected value is computed independently.
- Where a routine merely delegates, assert the *contract* — the returned arrays equal what the routine stored — rather than re-deriving numbers a different routine's own test already covers.

### Test code layout

Test code is organized by role, one file each:

| File | Role |
|---|---|
| `test_reference.f90` | Tolerances, reference values, and `ref_*` reimplementations computed independently of the routine under test. |
| `opentrustregion_unit_tests.f90` | The Fortran-level `test_*` functions for the core solver. |
| `opentrustregion_mock.f90` | Mocks of the module's own production routines that are bridged to the C interface (`solver`, `stability_check`), so a test can verify a routine invoked them correctly without running the real logic. |
| `c_interface_unit_tests.f90` | Tests for the `bind(C)` wrapper layer itself; defines its own local `bind(C)` mocks rather than using `c_interface_mock.f90`. |
| `c_interface_mock.f90` | `bind(C)`-signature mocks used exclusively by the Python interface tests, dynamically loaded via `ctypes` from `libotrtestsuite`. |

### Shared architectural patterns

**Settings types.** `solver_settings_type` / `stability_settings_type` extend the abstract `settings_type`, with an `init` type-bound procedure and an `initialized` flag. The owning routine checks `settings%initialized` on entry and calls `init` if needed; C/Python wrappers pre-populate defaults via `*_init()`/`__init__` so the flag is `.true.` by the time Fortran sees it. Default values are defined once in a `default_*_settings` parameter and replicated in the C `*_init()` function and the Python `Structure` defaults — keep them in sync.

**The Python → C → Fortran call chain:** `solver(...)` in `python_interface.py` builds a `SolverSettings` ctypes struct, wraps Python callbacks with `CFUNCTYPE`, and calls into `libopentrustregion`'s C entry point (`src/c_interface.f90`), which stores the C function pointers in module-level `procedure(..)_before_wrapping` slots and invokes `standard_solver` (the real `solver`) with Fortran-shaped wrappers that call back through those pointers. **Consequence:** the C interface is not re-entrant across threads — module-level pointer state is shared. A signature change must keep all layers (Fortran abstract interface → C wrapper + abstract interface → C typedef → Python `CFUNCTYPE` + `Structure`) in lockstep.

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

```sh
# full suite (Python driver calling into libotrtestsuite: Fortran unit + system tests)
python3 -m pyopentrustregion.testsuite                  # from an installed/editable build
python3 pyopentrustregion/testsuite.py                  # from the source tree against ./build

# single test class or method (stdlib unittest)
python3 -m unittest pyopentrustregion.testsuite.SystemTests
python3 -m unittest pyopentrustregion.testsuite.OpenTrustRegionTests.test_solver
```

The Python suite runs Fortran- and C-side tests (via symbols dynamically loaded from `libotrtestsuite`) and pure-Python wrapper tests. System tests need `pyopentrustregion/test_data/*.bin`.

**Where to add coverage for a new setting.** As assertions inside the existing unit test for the routine the setting affects (`tests/opentrustregion_unit_tests.f90`), not as a new system test — system tests exercise real chemistry data end-to-end, and one per settings combination would explode combinatorially. When a setting can push a routine down a materially different code path (interior vs. boundary-crossing branch in a trust-region solver), exercise both — a check that only hits the "easy" branch (e.g. always starting near a minimum) can pass while a branch reached only near a saddle point stays untested.

### Registering a new Fortran unit test

A `logical(c_bool) function test_<routine>() bind(C)` is only reachable once wired into both:

1. The source file in the CMake test source list (`OPENTRUSTREGION_TESTS`).
2. The `fortran_tests["opentrustregion_tests"]` (or `"c_interface_tests"`) name list in `pyopentrustregion/testsuite.py`, alphabetically ordered, without the `test_` prefix — the `@add_tests` decorator on the corresponding `unittest.TestCase` turns each name into a `test_<name>` method.

A name in the list with no matching Fortran symbol fails when that test runs, with `AttributeError: dlsym(...): symbol not found`.

**Exactly one registered test per production routine** — never `test_<routine>_<case_a>`, `test_<routine>_<case_b>`, … . Cover multiple configurations as cases *inside* that one test: loop where setup can be shared, or call a per-case `check_*` helper that is deliberately not `bind(C)` and not registered. Keep the case name in the failure message (`"test_<routine> failed for <case>: ..."`) since the harness only reports the registered name.

A `check_*` (or other non-`bind(C)`) helper earns its existence by eliminating real duplication, never merely to keep one case's code out of the registered test's body — write one function parameterized over the case, or a loop with a shared body, and reserve separate helper routines for pieces reused *across* cases (a shared fixture, a shared assertion). Place any such helper near the top of the file, right after the mocks and before the first `test_*` function.

Because a routine's cases now live behind one entry, a silently-skipped case is invisible in suite output — have the driving test run every case unconditionally and `and` the results together, rather than returning early on the first failure.

**Verify a new test actually fails when the routine is broken.** Mutate the routine (flip a sign, swap an index, drop a term), rebuild, confirm the test fails, then restore. Mutating several routines at once and checking the failure set matches one-to-one is efficient for a whole suite. Some mutations are legitimately benign for a given input (scaling safeguards, guards on paths the test doesn't reach) — pick a different mutation rather than concluding the test is weak.

**Before finishing:** build with `-DCMAKE_BUILD_TYPE=Debug` (`-Wall -Wextra -fcheck=all`) and run the suite — it catches out-of-bounds accesses release silently tolerates.

### CMake options that matter

- `INTEGER_SIZE` (`4` or `8`): library integer width. Unset → CMake autodetects, trying 32-bit BLAS/LAPACK first. Output library named `libopentrustregion_32.*` / `_64.*`. `USE_ILP64` (auto-set when `INTEGER_SIZE=8`) switches Fortran `ip` and C `c_ip` to 64-bit and remaps BLAS/LAPACK symbols to `_64` variants when `check_fortran_function_exists` finds them.
- `BLAS_LIBRARIES` / `LAPACK_LIBRARIES`: must be set together with `INTEGER_SIZE` if overriding autodetection — one without the other is a fatal error.
- `OpenTrustRegion_BUILD_TESTING` (default `ON` when top-level): builds `libotrtestsuite`, which Python loads to drive the Fortran tests.
- `CONDA_BUILD=1` env var: `setup.py` skips the embedded CMake invocation (the conda recipe builds the C library separately).

### Preprocessing and integer kinds

All Fortran sources compile with `Fortran_PREPROCESS ON`. Integer kind selection and the BLAS/LAPACK 64-bit symbol remap (`ddot=ddot_64`, etc.) happen via `#ifdef USE_ILP64` and `add_compile_definitions` in CMake — never hardcode integer kinds.

BLAS/LAPACK are called through implicit interfaces (`external :: dgemm`), so gfortran infers each dummy argument's kind from the first call site and rejects a later call that disagrees. Both rules below are invisible in the default 32-bit build, where `ip` is `int32` and `1_ip == 1`:

- **Every integer argument to a BLAS/LAPACK routine must be `ip`-kinded** — `1_ip` not `1` for increments/leading dimensions, `size(x, kind=ip)` when passing on a `size()` result. Same for integers passed to project routines with an explicit `integer(ip)` dummy. Locals used as LAPACK arguments must be declared `integer(ip)`.
- **A newly used BLAS/LAPACK symbol must be added to the remap lists in `CMakeLists.txt`** (`add_compile_definitions("name=name_64")`, guarded by `BLAS_64`/`LAPACK_64`) — including test sources (`zheev` reaches the build only via `tests/opentrustregion_system_tests.f90`). A missing entry links the 32-bit-integer symbol from an ILP64 build, corrupting arguments at runtime rather than failing to build.

**Verifying an ILP64 build without an ILP64 BLAS.** Most dev machines only have 32-bit-integer BLAS, so `check_fortran_function_exists("sgemm_64")` fails and the remap is never exercised. Force it by pre-seeding the cache variables and checking the symbols the objects actually reference:

```sh
cc -shared -o /tmp/blas64stub.dylib /tmp/stubs.c   # one empty function per name_64_ symbol
cmake -S . -B /tmp/build_ilp64 -DINTEGER_SIZE=8 \
      -DBLAS_LIBRARIES=/tmp/blas64stub.dylib -DLAPACK_LIBRARIES=/tmp/blas64stub.dylib \
      -DBLAS_64=1 -DLAPACK_64=1 -DBUILD_SHARED_LIBS=ON
cmake --build /tmp/build_ilp64
find /tmp/build_ilp64 -name '*.o' | xargs nm -u | grep -E '_(d|z)[a-z]+_'   # none may lack _64
```

A clean link proves nothing is left unmapped; it does **not** verify numerical behaviour, which needs a real ILP64 BLAS.

### Known gotchas

- **A `size()`/literal-kind mistake in a BLAS call passes the default build and only breaks under `INTEGER_SIZE=8`**, as `Error: Type mismatch between actual argument at (1) and actual argument at (2) (INTEGER(4)/INTEGER(8))` pointing at two unrelated call sites of the same routine — the two that disagree, not the one that's wrong. Fix by making every integer argument `ip`-kinded, not by changing the named site.
- **`gfortran -fsyntax-only -I<moddir>` gives false confidence.** It checks only the pointed-to file against whatever `.mod` files already exist — it doesn't re-verify those `.mod`s. Editing a `type` whose fields are used across modules can leave a stale consumer `.mod` "passing" syntax-only checks while a real build breaks with `Fatal Error: Mismatch in components of derived type '...': expecting 'X', but got 'Y'`. Always confirm interface changes with `cmake --build`, not `-fsyntax-only`.
- **A `build/` directory created by `pip install` can't be rebuilt directly later.** pip's ephemeral `cmake` path gets baked into `build/CMakeCache.txt` (`CMAKE_COMMAND`) and generated Makefile stamp rules. Once pip's temp env is gone, `cmake --build build` fails with `<temp-path>/cmake: No such file or directory`. Diagnose with `grep CMAKE_COMMAND build/CMakeCache.txt`. Fix: re-run `pip install -e .`, or maintain a separate manually-configured build directory with the system `cmake`.
- **The Python driver only looks for `../build`.** `python_interface.py` searches site-packages, then `pyopentrustregion/`, then `<repo>/../build` — a manually-configured `build_manual/` is invisible to it, and it silently loads whatever's in `build/` instead. To drive a custom build directory, run from a scratch directory with symlinks named `pyopentrustregion` and `build`:

  ```sh
  mkdir -p /tmp/run && cd /tmp/run
  ln -sfn <repo>/pyopentrustregion pyopentrustregion
  ln -sfn <repo>/build_manual build
  python3 -m pyopentrustregion.testsuite
  ```
- **`ref_settings_type_c` in `test_reference.f90` has an implicit, unenforced field-order convention.** `PyInterfaceTests.setUpClass` in `pyopentrustregion/testsuite.py` rebuilds a matching struct by walking `SolverSettings.c_struct._fields_ + StabilitySettings.c_struct._fields_` (that concatenation order, bools-then-reals-then-ints, then character arrays in relative order) and passes it into Fortran's `get_reference_values` via a raw `POINTER(RefSettingsC)` cast. So `ref_settings_type_c` must declare all `solver_settings_type_c`-owned character/string fields before all `stability_settings_type_c`-owned ones, regardless of what feels natural. Getting it wrong causes no crash — the two structs read across each other's byte ranges, so a field silently reports another field's value, or raises `UnicodeDecodeError` if the misaligned bytes aren't valid UTF-8. A new settings field failing with an inexplicable value mismatch or decode error → suspect field order here first.

## Core library

### Source layout

- `src/opentrustregion.f90` — the numerical core: `solver`, `stability_check`, the Davidson/Jacobi-Davidson/truncated-CG subsystem solvers, settings derived types, callback abstract interfaces, error-code constants. Single module, several thousand lines.
- `src/c_interface.f90` — `bind(C)` wrapper module. Stores Fortran procedure pointers to C callbacks at module scope (`update_orbs_before_wrapping`, etc.) and adapts C-style `(*)` arrays + return-code functions into Fortran-style `(:)` arrays + `intent(out) :: error` subroutines.
- `include/opentrustregion.h` — C header mirroring `solver_settings_type` / `stability_settings_type` as C structs, plus `*_init()` helpers and `solver`/`stability_check` prototypes.
- `pyopentrustregion/python_interface.py` — ctypes wrapper. `SolverSettings`/`StabilitySettings` as `ctypes.Structure` mirrors of the C structs, Python callbacks wrapped with `CFUNCTYPE`, integer error codes converted to `RuntimeError`.
- `tests/` — `opentrustregion_unit_tests.f90` / `c_interface_unit_tests.f90` (Fortran/C interface unit tests, both in `libotrtestsuite`), `opentrustregion_system_tests.f90` (system tests against `pyopentrustregion/test_data/`), `c_system_tests.c` (system tests through the C interface), `*_mock.f90` (mock callbacks for both unit suites), `test_reference.f90` (tolerances + reference values).

### Fortran/C/Python interfaces must stay consistent

The same callback signatures, settings fields, and defaults are described in seven places that must agree — any change to a callback signature, a settings field (add/remove/reorder), or a default value must land in all seven in the same PR:

1. Fortran abstract interfaces and `solver_settings_type`/`stability_settings_type` in `src/opentrustregion.f90`.
2. C abstract interfaces and `bind(C)` `solver_settings_type_c`/`stability_settings_type_c` in `src/c_interface.f90`.
3. C struct layouts/typedefs in `include/opentrustregion.h`.
4. `SolverSettingsC`/`StabilitySettingsC` ctypes `_fields_` and `CFUNCTYPE` declarations in `pyopentrustregion/python_interface.py`.
5. `default_solver_settings`/`default_stability_settings` (Fortran) ↔ `solver_settings_init`/`stability_settings_init` (C) ↔ the Python `Settings` wrapper defaults.
6. Argument lists and snippets in `README.md`.
7. `ref_settings_type`/`ref_settings_type_c` in `tests/test_reference.f90`, used by the settings round-trip tests in `pyopentrustregion/testsuite.py` — see the field-*order* gotcha above.

Error-origin codes (`error_obj_func` etc.) must also stay synchronized with the README error-code table.

### Error codes

Encoded as `OOEE` integers (origin × 100 + specific code, see README). The Fortran core uses `add_error_origin` to tag a non-zero callback error with its origin (e.g. `error_update_orbs = 1200`). C entry points return the integer directly; Python wrappers raise a `RuntimeError` that includes the raw code as-is — they don't decode it into origin/specific parts, so callers must read the README table themselves. New origins: add as a parameter in `opentrustregion.f90` and keep the README table in sync.
