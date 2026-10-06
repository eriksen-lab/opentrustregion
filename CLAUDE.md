# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

OpenTrustRegion is a Fortran library implementing a second-order trust-region optimizer, exposed via Fortran, C, and Python (ctypes) interfaces — the same compiled shared/static library is consumed by all three.

## General conventions

**Formatting.** Fortran lines max 88 columns. Layout, spacing, keyword spelling, string literals, comments and blank lines of Fortran code are defined by `tools/f90_layout.py`, whose docstring is the authoritative rule list; don't restate its rules here or apply them by hand. Run `python3 tools/f90_layout.py` after every Fortran edit: by default it checks only statements touched relative to `HEAD`, `--fix` applies the fixes, and `--all` checks whole files (with `--fix` it requires explicit paths, e.g. `git ls-files '*.f90' | xargs python3 tools/f90_layout.py --all --fix`). CI runs `--all` on every tracked `.f90` file. The default mode also checks every untracked, non-ignored `.f90` file in full, so a bare `--fix` rewrites scratch copies of Fortran files left in the tree too; move them out first. Statements containing a comment, a `;` or a continuation line starting with `&` are never reflowed (their spacing and indentation are still fixed), so shorten an overlong one by hand. Every procedure opens with a `!`-delimited comment block describing what it does; each logical step inside gets a lowercase `!` comment. Unit tests: assume success (`test_<name> = .true.`), then one `if (...) then` / `write(stderr, *) "test_<name> failed: ..."` / `test_<name> = .false.` block per assertion, no early returns except where a later assertion would crash.

**Python formatting.** All Python (`pyopentrustregion/`, `setup.py`, `tools/`) is `black`-formatted, default settings. Run `black pyopentrustregion setup.py tools` before considering Python changes done.

**C formatting.** All C headers and test sources (`include/`, `tests/*.c`) are `clang-format`-formatted per the repo-root `.clang-format` (LLVM style, 88-column limit to match the Fortran convention above). Run `clang-format -i` on touched C/H files before considering C changes done.

**Clarity over performance.** The library is not on the hot path of a quantum-chemistry calculation — the host program's integral transforms and Hessian-vector products dominate. Prefer short, obviously-correct code over fast code. Don't propose performance refactors (buffer growth, pooling, micro-optimizations) without evidence the affected code is hot for a real workload.

**Norms: BLAS in production, intrinsics in tests.** Production Fortran (`src/`) uses BLAS (`dnrm2`, `ddot`) for norms/dot products. Test code (`tests/`) uses intrinsics (`norm2`, `dot_product`, `matmul`, `sum`, `transpose`) instead. BLAS/LAPACK is allowed in tests only where no intrinsic exists (`dsyev`, `dgeev`, `zheev`) — don't hand-roll linear algebra to avoid the dependency.

### Unit tests must not depend on untested code

A unit test may only call the routine it is testing — never another production routine to build its inputs or expected values, or the test silently inherits that routine's correctness.

- Construct inputs directly in the test, randomly or with LAPACK, not by calling the production routine that would normally produce them.
- When a routine under test internally calls another module's routines, reimplement those as `ref_*` helpers so the expected value is computed independently.
- Where a routine merely delegates, assert the *contract* — the returned arrays equal what the routine stored — rather than re-deriving numbers a different routine's own test already covers.
- Exercise only the options the routine itself distinguishes. Where it calls another production routine, its test checks only that the callee was called (with the right arguments, and that its result is used correctly), not the callee's options or the input regimes that drive the callee down different branches; those belong in the callee's own test. E.g. `test_solver` checks that the solver, with `line_search` set, hands the step to `bracket`, while which bracketing branch an objective function takes is covered in `test_bracket`.

### Test code layout

Test code is organized by role, one file each:

| File | Role |
|---|---|
| `test_reference.f90` | Tolerances, reference values, and `ref_*` reimplementations computed independently of the routine under test. The reference settings `ref_solver_settings`/`ref_stability_settings` are values of the production types with a distinct value in every field but the logicals, which the tests set one at a time instead; C and Python get them field by field without the conversion routines under test, as a filled C struct (`get_reference_solver_values`) or by field name (`reference_field`), so their tests hold no values of their own. `callbacks_unset`, `unset_callbacks` and `callbacks_wrapped` check, unset or check the wrapping of every callback function of Fortran or C settings, so a new callback function is listed there once instead of in every test. The checkers `check_<callback>_funptr`/`check_<callback>_c_funptr`, which the mocks use to check the function pointers they receive, sit side by side as Fortran/C pairs. Also the host-context fixtures shared by both unit suites: `host_context_type` with its `host_context` and `stability_host_context` instances, `arm_host_context`/`arm_host_context_c`, `check_host_context`/`check_host_context_c` and `host_context_reached`. |
| `opentrustregion_unit_tests.f90` | The Fortran-level `test_*` functions, plus mocks/fixtures needed only to drive them. |
| `opentrustregion_mock.f90` | Mocks of the module's own production routines that are bridged to the C interface (`solver`, `stability_check`), so a test can verify a routine invoked them correctly without running the real logic. |
| `c_interface_unit_tests.f90` | Tests for the `bind(C)` wrapper layer itself; defines its own local `bind(C)` mocks rather than using `c_interface_mock.f90`. |
| `c_interface_mock.f90` | `bind(C)`-signature mocks used exclusively by the Python interface unit tests, dynamically loaded via `ctypes` from `libotrtestsuite`. |

**`ref_*` naming.** The prefix marks a role — an independent reimplementation of what a production routine computes — not a location; it travels with the function wherever the placement rule moves it, and it is decided by what the helper computes, not by whether a call site compares against it or feeds it in as an input. A fixture that builds a *plausible* input (randomly, via LAPACK, or by reassembling a routine's own output) does not take the prefix, even in the same file. Either kind may keep an atypical name purely to dodge a collision with a production routine or local variable in scope at its call site — check for a collision before insisting on the "cleaner" name.

**Module-level `use` in test files is restricted.** In `<name>_unit_tests.f90` and siblings, module-level `use` statements above `contains` are limited to `rp`, `ip`, their C counterparts `c_rp` and `c_ip`, `kw_len`, `stderr`, `stdout`, `tol`, `tol_c`, the dimension `n_param`/`n_param_c` and the reference settings `ref_solver_settings`/`ref_stability_settings` (with the operators comparing against them) that `test_reference.f90` fixes for the C-interface fixtures, and intrinsic module bindings (`iso_c_binding` and the like) — check other modules for the exact set before adding to it. Names from `iso_c_binding` are always imported at module level, never inside a procedure. The exceptions are what a declaration in the module's specification part needs: the parent of a test-local type extension, since a type with type-bound procedures can only be declared there, the interface of a procedure pointer that checks a mock complies with it (`procedure(update_orbs_c_type), pointer :: mock_update_orbs_ptr => mock_update_orbs`), and the type of a variable in which a mock records what it received (`received_stability_callbacks` in `opentrustregion_mock.f90`). Everything else (mock callbacks, `setup_settings`, other shared dimension parameters, `ref_*` helpers) is imported locally inside the specific `test_*` function that needs it, even if several functions in the file need the same symbol. A module-level `use` beyond this set pollutes every procedure in the file and hides which test actually relies on what.

**Prefer shared dimension-parameter fixtures over local literals**, unless it complicates the test a lot. The test reference modules define canonical dimensions that most unit tests should import rather than redeclaring. The core unit tests' `n_param` is instead the dimension of the Hartmann 6D problem most of them share, defined with that problem in `opentrustregion_unit_tests.f90`; `test_reference.f90`'s `n_param` is the dimension of the C-interface fixtures.
- A test whose local variable already carries the name renames the import on its `use` line only, with the suffix `_ref`, never the prefix `ref_`, which marks independent reimplementations.
- Exception: a test whose expected values are hand-derived for a specific small/structured case (a closed-form 2×2 rotation, 3×3 matrices with hand-computed results, a matrix crafted to trigger a rank-detection path) may keep local literals and declare its own dimensions, since forcing the shared size would mean re-deriving the numbers. When correctness is checked via algebraic properties rather than exact values, use the shared fixture.

**Prefer generated data over hand-typed literals** (`generate_random_symm_matrix` / `call random_number`), unless the specific values are load-bearing. A hand-typed `reshape([...])` is justified only when the test checks a specific hand-derivable closed-form target depending on those exact numbers. When every assertion is relational (a copy, a fixed scaling, an algebraic identity, or a value the test recomputes independently from the same input) the specific values never mattered, and random generation exercises the same code path while being harder to accidentally pass with a bug. Before hand-typing a matrix meant to be "just some valid state" (often signaled by the test calling it "arbitrary"), check whether random generation works instead. Where a matrix must stay block-diagonal/diagonal/structurally constrained to match what production produces, keep that structure but randomize within it.

### Shared architectural patterns

**Settings types.** `solver_settings_type` / `stability_settings_type` extend the abstract `settings_type`, with an `init` type-bound procedure and an `initialized` flag. The owning routine checks `settings%initialized` on entry and calls `init` if needed; `solver` does the same for its nested `stability_settings` up front, since the internal stability check's own `init` would discard the callbacks and context the solver hands down to it. C/Python wrappers pre-populate defaults via `*_init()`/`__init__` so the flag is `.true.` by the time Fortran sees it. C settings that arrive uninitialized anyway are replaced by the defaults in the `assign_*_f_c` conversion rather than left to `init`, which would also discard the callback bundle the C wrapper stores in `settings%context`; as in Fortran, such settings carry no host context, so `store_optional_c_callbacks` skips the optional callbacks and context of uninitialized C settings. Default values are defined once in a `default_*_settings` parameter; the C `*_init()` helpers and the Python `Settings` wrappers obtain them by calling `init_*_settings_c` (C name `init_*_settings`) in `src/c_interface.f90`.

**The host context.** Every callback interface takes `class(*), intent(in), pointer :: context` as its last dummy (`void *context` in C), and `settings_type` carries a matching `class(*), pointer :: context` component. The library passes `settings%context` at every call site unconditionally — a *pointer* dummy accepts a disassociated actual, which is why no call site branches on `associated`, and why the component must be `pointer` rather than `allocatable` (the internals declare `settings` `intent(in)`, which protects the pointer association but not the target). `solver` (in `check_stationary_point`) lends its context to `settings%stability_settings%context` when that is unset, for the duration of the internal stability check only, and takes it back afterwards. When the nested context *is* set it wins: the internal check hands it to every callback it calls, including the `hess_x` returned by the solver's `update_orbs` and the inherited `precond`/`project`/`logger`, so a host that sets one must make those callbacks accept it (documented in the README). `test_check_stationary_point` checks that the nested settings keep their context alongside the inherited callbacks and hand it to `hess_x`, `test_stability_check` that every callback receives `settings%context`, and `test_solver_c_wrapper` that the C wrapper's nested bundle holds the nested settings' own callbacks and context where they were provided and initialized and the solver's otherwise. Unlike `precond`/`project`/`logger`, which stay inherited on the host's settings object, a data pointer left there would carry over to the next solve on a reused settings object and could dangle; `test_check_stationary_point` asserts nothing is left behind.

**The Python → C → Fortran call chain:** `solver(...)` in `python_interface.py` takes the `SolverSettingsC` ctypes struct a `SolverSettings` object wraps, wraps Python callbacks with `CFUNCTYPE`, and calls into `libopentrustregion`'s C entry point (`src/c_interface.f90`), which collects the C function pointers plus the host's `void*` in a `c_callbacks_type` bundle on the entry point's own stack frame, stores that bundle in `settings%context`, and invokes `standard_solver` (the real `solver`) with Fortran-shaped wrappers that `select type` the bundle back out. **Consequence:** there is no module-level callback state, so solves nest and run concurrently. The one piece of process-wide state left is the Fortran runtime's random number generator, which `init_rng` reseeds on every `solver`/`stability_check` call, so a nested or concurrent call perturbs the other's random trial vectors, not its correctness. The internal stability check gets a *second* bundle that starts as a copy of the solver's (so inherited callbacks resolve) and is then overridden by whatever `settings_c%stability_settings` provides — a nested bundle whose slot is empty while the Fortran pointer is inherited would call a null procedure pointer. A signature change must keep all layers (Fortran abstract interface → C wrapper + abstract interface → C typedef → Python `CFUNCTYPE` + `Structure`) in lockstep. Python needs no context: `python_interface.py` accepts arbitrary callables, so closures already carry host data.

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
# full suite (Python driver calling into libotrtestsuite: unit, integration and system tests)
python3 -m pyopentrustregion.testsuite                  # from an installed/editable build
python3 pyopentrustregion/testsuite.py                  # from the source tree against ./build

# single test class or method (stdlib unittest)
python3 -m unittest pyopentrustregion.testsuite.OpenTrustRegionSystemTests
python3 -m unittest pyopentrustregion.testsuite.OpenTrustRegionUnitTests.test_solver
```

The Python suite runs Fortran- and C-side tests (via symbols dynamically loaded from `libotrtestsuite`) and pure-Python wrapper tests. System tests need `pyopentrustregion/test_data/*.bin`.

**Where to add coverage for a new setting.** As assertions inside the existing unit test for the routine the setting affects (`tests/opentrustregion_unit_tests.f90`), not as a new system test — system tests exercise real chemistry data end-to-end, and one per settings combination would explode combinatorially. When a setting can push a routine down a materially different code path (interior vs. boundary-crossing branch in a trust-region solver), exercise both — a check that only hits the "easy" branch (e.g. always starting near a minimum) can pass while a branch reached only near a saddle point stays untested.

**What the system and integration tests cover.** The tests come in three tiers: unit tests run one routine against mocks, integration tests run the real library through one interface's declarations, and system tests run the real library on real data. Each tier checks what it adds, not the scenarios of the tier below. The system tests check numerical behaviour:

| Test | Covers |
|---|---|
| `test_h2o_fb_solver_<option>` (`opentrustregion_system_tests.f90`) | `solver` on the water Foster-Boys localization with the default settings and with one algorithm-selecting option switched on at a time (Jacobi-Davidson, TCG, `line_search`, `stability`), each from a guess and from a saddle point; checks the converged gradient, the reference minimum and its stability. |
| `test_h2o_fb_stability_check_<option>` (same file) | `stability_check` with the default Davidson and with the Jacobi-Davidson diagonalization solver at the minimum (stable) and at a saddle point (unstable, returned direction along the eigenvector of the lowest Hessian eigenvalue, since the saddle point has several unstable modes). |

The integration tests check what crosses the C and Python interfaces. The interface unit tests replace the other side with mocks, so these are the only tests that run the real entry points behind the real header and ctypes declarations; their solves on the Hartmann 6D problem only serve as evidence that the results came back:

| Test | Covers |
|---|---|
| `test_solver_c`, `test_stability_check_c` (`c_integration_tests.c`) | The real library through the C header: every optional callback of both settings structs, including the nested stability settings' own, the host context reaching every callback (including those of the internal stability check, whose nested bundle carries the host's pointer), a stability check nested in a running solve, the returned counters and stability flag, and the call without a direction. |
| `test_settings_layout` (same file, and `PyIntegrationTests` for the ctypes structures in `python_interface.py`) | That every field of the C settings structs, including the nested stability settings, is read back under its own name from settings whose fields Fortran sets one by one to the reference values (`get_reference_solver_values` in `test_reference.f90`), compared against the value `reference_field` returns for the field's name, so any drift between the header and the `bind(C)` types is caught. The Python test loops over the ctypes fields, so a new field needs no change there. |
| `test_solver_settings_init`, `test_stability_settings_init` (same file) | That the header's inline `*_init()` wrappers return the library's defaults with no callbacks and no host context. Both sides come from the same Fortran assignment, so these do not detect layout drift; `test_settings_layout` does. |
| `test_solver_py`, `test_stability_check_py` (`PyIntegrationTests` in `testsuite.py`) | The real library through the Python wrapper: every optional callback of both settings classes, including the nested stability settings' own, the returned counters and stability flag, and the call without a direction. |

System tests are registered per option, not per routine: the one-registered-test-per-routine rule and the advice above to cover a new setting in the unit tests apply to unit tests and to settings that only tune a routine. A new option that selects a different algorithm gets its own registered `test_h2o_fb_solver_<option>` or `test_h2o_fb_stability_check_<option>` test, never a combination with another option. The integration tests gain a check only when a new callback, settings field or return value crosses their interface.

### Registering a new Fortran unit test

A `logical(c_bool) function test_<routine>() bind(C)` is only reachable once wired into both:

1. The source file in the CMake test source list (`OPENTRUSTREGION_TESTS`).
2. The `fortran_tests["opentrustregion_unit_tests"]` (or `"c_interface_unit_tests"`) name list in `pyopentrustregion/testsuite.py`, alphabetically ordered, without the `test_` prefix — the `@add_tests` decorator on the corresponding `unittest.TestCase` turns each name into a `test_<name>` method.

A name in the Python list with no matching Fortran symbol fails at import with `AttributeError: dlsym(...): symbol not found`.

System tests and C integration tests are wired the same way: their source files are in the same CMake list, and their names go into `fortran_tests["opentrustregion_system_tests"]` and `fortran_tests["c_integration_tests"]`.

**Exactly one registered test per production routine** — never `test_<routine>_<case_a>`, `test_<routine>_<case_b>`. Cover multiple configurations as cases *inside* that one test: loop where setup can be shared, or call a per-case `check_*` helper that is deliberately not `bind(C)` and not registered. Keep the case name in the failure message (`"test_<routine> failed: ... for <case>."`) since the harness only reports the registered name.

Because a routine's cases live behind one entry, a silently-skipped case is invisible in suite output — have the driving test run every case unconditionally and `and` the results together, rather than returning early on the first failure.

**Check that an error message was printed, not its text.** When the routine under test detects an error itself, assert the error code and, separately, that an error message was printed, with only error messages passed to the logger and the log cleared beforehand (`setup_error_logging` in `opentrustregion_unit_tests.f90`). Since the routine returns right after printing it, a non-empty log can only be that message; don't compare its text. An error a routine merely propagates from a failing callback prints nothing, so there is no message to check. Warnings keep their exact-text comparison: the routine continues after a warning, so other messages can follow it in the same log. Compare against the `*_warning_msg` parameter the production module defines for each warning instead of repeating its text in the test.

**Error propagation is covered by fault injection.** A routine's own test checks the origin it tags on an error at its own callback call sites (`error_hess_x + 1` for a failing Hessian linear transformation it calls), not the `if (error /= 0) return` lines that only pass on an error of an internal callee. Those are covered by the fault-injection sweeps in `test_solver` and `test_stability_check`: an unfaulted run counts every callback call in a `hartmann6d_fault_context_type`, then every one of these calls fails in turn and the test checks that the error is reported with the failing callback's origin and that the reported counters still agree with the calls. A new callback gets a `faulty_*` mock and an entry in the fault lists, and a new option that changes which callbacks the solver or stability check calls gets a case in their sweep.

**Helper procedures in test files.**
- A `check_*` (or other non-`bind(C)`) helper earns its existence by eliminating real duplication, never merely to keep one case's code out of the registered test's body — five near-identical one-case functions just relocate the duplication; write one function parameterized over the case, or a loop with a shared body.
- A helper that would be called once without arguments, only reading and writing its host's variables, is part of the host's body and stays inline. A single call justifies a separate procedure only for a mock callback, which has to be a procedure, or for a function of its arguments used in an expression, such as a `ref_*` reimplementation.
- A helper only one procedure uses is an internal procedure of that procedure, after its own `contains`, unless it is one of a family of `ref_*` twins whose other members live at module level. An internal procedure cannot contain another, so helpers calling each other are nested side by side in their common caller. Mocks in `opentrustregion_mock.f90`/`c_interface_mock.f90` stay module procedures however many callers they have, since those files exist to export them.
- Mock callbacks, which a test assigns to a procedure pointer, likewise stay module procedures however few callers they have: gfortran implements a pointer to an internal procedure with a trampoline, which on Linux needs an executable stack that glibc ≥ 2.41 refuses to `dlopen`, breaking the ctypes-driven testsuite (macOS uses heap trampolines, so it looks fine locally).
- An internal procedure uses its host's variables and imports directly instead of redeclaring them or taking them as arguments that every call fills with the host's same-named variable. What must stay its own (loop indices, which would otherwise overwrite a host loop calling it, and dummies receiving a slice or a different array) gets a name the host does not use, so that it never masks a host name.
- Helpers used by several procedures sit together near the top of the file, before the first `test_*` function, grouped by role so that each group only uses the ones above it: mocks and their support, then generators of random or structured data, then fixtures that set up objects or collect their state, then `ref_*` reimplementations, then `check_*` helpers. Where a group uses helpers of a group that this order puts after it, as the reference-settings fixtures in `test_reference.f90` use `ref_character_to_c`, the dependency order wins and the used group moves up. Within a group, helpers concerning the same object stay next to each other, in the order of the production code they serve (deterministic constructors ahead of random generators); `check_*` helpers follow the order of the tests they drive.
- Names follow the role: `generate_random_<what>` for a generator of random data, `<what>_matrix` for a deterministic constructor (`identity_matrix`), `setup_<what>` for a fixture that sets up an object.

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

## Core library

### Source layout

- `src/opentrustregion.f90` — the numerical core: `solver`, `stability_check`, the Davidson/Jacobi-Davidson/truncated-CG subsystem solvers, settings derived types, callback abstract interfaces, error-code constants. Single module, several thousand lines.
- `src/c_interface.f90` — `bind(C)` wrapper module. Collects the C callback pointers in a per-call `c_callbacks_type` bundle held on the entry point's stack and passed through as the opaque context, and adapts C-style `(*)` arrays + return-code functions into Fortran-style `(:)` arrays + `intent(out) :: error` subroutines. The module-scope `solver`/`stability_check` procedure pointers remain as test-injection seams; only the wrappers' unit tests point them at mocks, restoring them afterwards, and nothing sets them per call.
- `include/opentrustregion.h` — C header mirroring `solver_settings_type` / `stability_settings_type` as C structs, plus `*_init()` helpers and `solver`/`stability_check` prototypes.
- `pyopentrustregion/python_interface.py` — ctypes wrapper. `SolverSettings`/`StabilitySettings` wrapper classes around the `SolverSettingsC`/`StabilitySettingsC` `ctypes.Structure` mirrors of the C structs, Python callbacks wrapped with `CFUNCTYPE`, integer error codes converted to `RuntimeError`.
- `tests/` — `opentrustregion_unit_tests.f90` / `c_interface_unit_tests.f90` (Fortran/C interface unit tests, both in `libotrtestsuite`), `opentrustregion_system_tests.f90` (system tests against `pyopentrustregion/test_data/`), `c_integration_tests.c` (integration tests of the real library through the C header), `opentrustregion_mock.f90` (mocks of `solver`/`stability_check` for the C-interface unit tests), `c_interface_mock.f90` (`bind(C)` mocks for the Python interface unit tests), `test_reference.f90` (tolerances, reference values and shared fixtures). The C and Python integration tests read the Hartmann 6D problem from the `hartmann6d_*` `bind(C)` data that `opentrustregion_unit_tests.f90` exports.

### Fortran/C/Python interfaces must stay consistent

The same callback signatures, settings fields, and defaults are described in seven places that must agree — any change to a callback signature or a settings field (add/remove/reorder) must land in all seven in the same PR, while a default value changes only in its definition (item 5) and has to stay distinct from its reference value (item 7):

1. Fortran abstract interfaces and `solver_settings_type`/`stability_settings_type` in `src/opentrustregion.f90`.
2. C abstract interfaces and `bind(C)` `solver_settings_type_c`/`stability_settings_type_c` in `src/c_interface.f90`.
3. C struct layouts/typedefs in `include/opentrustregion.h`.
4. `SolverSettingsC`/`StabilitySettingsC` ctypes `_fields_` and `CFUNCTYPE` declarations in `pyopentrustregion/python_interface.py`.
5. The default values in `default_solver_settings`/`default_stability_settings` in `src/opentrustregion.f90`, their only definition: the C `solver_settings_init`/`stability_settings_init` helpers and the Python `Settings` wrappers fill their structs by calling `init_solver_settings_c`/`init_stability_settings_c` (C names `init_solver_settings`/`init_stability_settings`) in `src/c_interface.f90`.
6. Argument lists and snippets in `README.md`.
7. The reference settings in `tests/test_reference.f90`: the value in `ref_solver_settings`/`ref_stability_settings`, the comparisons `equal_*` and `equal_*_c`, the field-by-field export `get_reference_*_values` and the lookup `reference_field`, plus its check in the C `test_settings_layout` (a `check_field`/`check_keyword` line, or for a logical its name and field in the two parallel arrays a static assert keeps in step), and for a callback its line in `callbacks_unset`, `unset_callbacks` and `callbacks_wrapped`, a distinct address in `get_reference_*_values` and its `_Static_assert` offset at the top of `c_integration_tests.c`. The Python settings and layout tests read the fields by name and need no change. A new or changed callback also reaches the tests outside these seven places: its checkers `check_<callback>_funptr` and `check_<callback>_c_funptr` in `test_reference.f90`, the mocks of both unit suites including a `faulty_*` one for the fault-injection sweeps, the callbacks of `c_integration_tests.c` and the callback names listed in the Python tests in `testsuite.py`.

Error-origin codes (`error_obj_func` etc.) must also stay synchronized with the README error-code table.

### Error codes

Encoded as `OOEE` integers (origin × 100 + specific code, see README). The Fortran core uses `add_error_origin` to tag a non-zero callback error with its origin (e.g. `error_update_orbs = 1200`). C entry points return the integer directly; Python wrappers raise a `RuntimeError` that includes the raw code as-is — they don't decode it into origin/specific parts, so callers must read the README table themselves. New origins: add as a parameter in `opentrustregion.f90` and keep the README table in sync.
