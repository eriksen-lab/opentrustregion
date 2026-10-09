![License](https://img.shields.io/github/license/eriksen-lab/opentrustregion)
![CI](https://github.com/eriksen-lab/opentrustregion/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/github/eriksen-lab/opentrustregion/graph/badge.svg?token=NJIM6FDADD)](https://codecov.io/github/eriksen-lab/opentrustregion)
[![DOI](https://zenodo.org/badge/855894860.svg)](https://doi.org/10.5281/zenodo.17142554)

# OpenTrustRegion: A Reusable Library for Second-Order Trust Region Orbital Optimization

This library provides a robust and flexible implementation for second-order trust region orbital optimization, with extensive customization options to suit various use cases.

The following paper documents the theory and implementation of the methodology in OpenTrustRegion, and should be cited in any work using OpenTrustRegion:  

- Greiner, J.; Høyvik, I.-M.; Lehtola, S.; Eriksen, J. J. 
  A Reusable Library for Second-Order Orbital Optimization Using the Trust Region Method. Journal of Chemical Theory and Computation 2026, 22(2), 881–895. 
  DOI: [10.1021/acs.jctc.5c01576](https://doi.org/10.1021/acs.jctc.5c01576). 
  arXiv: [2509.13931](https://arxiv.org/abs/2509.13931).

## Installation

### Fortran or C Usage
To build the library for Fortran or C:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The installation can be tested by running the ```testsuite.py``` file in the ```pyopentrustregion``` directory.

### Python Usage
To install the library for use with Python:

```sh
pip install .
```

The installation can be tested by running

```sh
python3 -m pyopentrustregion.testsuite
```

### CMake Configuration Options

The build process can be customized using the following CMake options:

| Option | Type | Default | Description |
|--------|------|----------|------------|
| **BUILD_SHARED_LIBS** | `BOOL` | `OFF` | Build shared libraries (`.so`, `.dylib`) instead of static ones. |
| **OpenTrustRegion_BUILD_TESTING** | `BOOL` | `ON` | Build the project’s testsuite. |
| **OpenTrustRegion_INSTALL_CMAKEDIR** | `STRING` | (auto) | Project install directory. |
| **CMAKE_BUILD_TYPE** | `STRING` | `Release` | Choose the build type (`Debug`, `Release`, etc.). |
| **INTEGER_SIZE** | `STRING` | *(auto)* | Set the integer precision to `4` (32-bit) or `8` (64-bit). Required when providing custom BLAS/LAPACK libraries. Otherwise defaults to 32-bit integers and tries to locate compatible BLAS and LAPACK libraries. Falls back to 64-bit integers if 32-bit libraries cannot be found. The resulting library name reflects the chosen integer precision (`libopentrustregion_32.*` or `libopentrustregion_64.*`) |
| **BLAS_LIBRARIES** | `PATH` | *(auto)* | Path(s) to BLAS libraries. If not provided, CMake attempts to locate a suitable BLAS automatically. |
| **LAPACK_LIBRARIES** | `PATH` | *(auto)* | Path(s) to LAPACK libraries. If not provided, CMake attempts to locate a suitable LAPACK automatically. |
| **OpenTrustRegion_HOST_PROVIDES_BLAS** | `BOOL` | `OFF` | When enabled, OpenTrustRegion will not attempt to detect or link BLAS/LAPACK and the testsuite is automatically disabled. The calling program must provide BLAS/LAPACK routines that expose the unsuffixed symbol names (for example `ddot`, `dsyev`) with an integer width matching `INTEGER_SIZE`; otherwise linking will fail loudly. |
| **OpenTrustRegion_ENABLE_XHOST** | `BOOL` | `ON` | Optimize the Release build for the current machine's instruction set (`-march=native` for GNU, `-xHost` for Intel). Automatically disabled when cross-compiling, regardless of this setting. Turn off when building for a different machine than the one compiling (e.g. packaging/conda builds). |

### Nested and Concurrent Calls

The library keeps no state of its own between calls, so a solver or stability check call can be nested inside a callback function of another call, or run concurrently with another call in a different thread, as long as each call is given its own settings object. A nested call re-enters procedures of the library that are still running, which Fortran only permits for procedures compiled for recursion, and a concurrent call needs their local variables on its own thread's stack. The library then has to be built with `-frecursive` for gfortran (implied by `-fopenmp`) or `-recursive` for Intel:

```sh
cmake .. -DCMAKE_Fortran_FLAGS=-frecursive                  # Fortran or C
CMAKE_FLAGS="-DCMAKE_Fortran_FLAGS=-frecursive" pip install .  # Python
```

Every call also reseeds the random number generator of the Fortran runtime, which the whole program shares, so a nested or concurrent call changes the random trial vectors, and with them the iterations, of the other call, though not its correctness.

## Program Interfaces

The PySCF interface is available as an extension hosted at https://github.com/eriksen-lab/pyscf_opentrustregion. To install it, simply add its path to the **`PYSCF_EXT_PATH`** environment variable:
```sh
export PYSCF_EXT_PATH=path/to/pyscf_opentrustregion
```
Usage examples can be found in the **`examples`** directory of the PySCF interface repository.

The interface supports Hartree–Fock and DFT calculations via the **`mf_to_otr`** function, which wraps PySCF **`HF`** and **`KS`** objects into their OpenTrustRegion counterparts. Similarly, localization methods are available through the **`BoysOTR`**, **`PipekMezeyOTR`**, and **`EdmistonRuedenbergOTR`** classes, and state-specific CASSCF calculations are supported via the **`casscf_to_otr`** function applied to a PySCF **`CASSCF`** object. All returned objects are fully compatible with the original PySCF classes and can be used interchangeably.

Optional settings can be adjusted by modifying object attributes directly. Orbital optimization and internal stability analysis are performed using the **`kernel`** and **`stability_check`** member functions, respectively.

## Usage

The optimization process is initiated by calling a `solver` subroutine. This routine requires the following input arguments:

### Required Arguments

- **`update_orbs`** (subroutine):  
  Accepts and applies a variable update (e.g., orbital rotation), updates the internal state, and provides:
  - Objective function value (real)
  - Gradient (real array, written in-place)
  - Hessian diagonal (real array, written in-place)
  - A **`hess_x`** subroutine that performs Hessian-vector products:
    - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
    - Returns an integer error code (0 for success, positive integers < 100 for errors)
    - Receives the host context as its last argument.
    - Has to be provided whenever the orbital update succeeds, otherwise the solver fails. In Fortran, the procedure pointer argument is therefore `intent(inout)` rather than `intent(out)`, since the solver disassociates it before every call to detect a missing one.
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`obj_func`** (function):  
  Accepts and applies a variable update (e.g., orbital rotation) and returns:
  - Objective function value (real)
  - An integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`settings`** (settings_type): Settings object which controls optional arguments as described below.

---

The following Fortran snippet demonstrates how to use the `solver` interface:

```fortran
use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type, solver_settings_type, solver

procedure(update_orbs_type), pointer :: update_orbs_funptr
procedure(obj_func_type), pointer :: obj_func_funptr
integer(ip) :: n_param, error
type(solver_settings_type) :: settings

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_macro = 100
settings%subsystem_solver = "tcg"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run solver
call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)

! read back output fields
print *, "Number of orbital updates:", settings%n_update_orbs
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized using the `init()` method of the derived type, and default settings can be overridden (here, `conv_tol` and `n_macro`).
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- Finally, the `solver` is called with the initialized settings and callback functions.
- After the call, output fields on `settings` (here, `n_update_orbs`) are populated and can be read like any other component.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdio.h>
#include <string.h>
#include "opentrustregion.h"

c_int n_param;

// set callback function pointers to existing implementations
update_orbs_fp update_orbs_funptr = update_orbs;
obj_func_fp obj_func_funptr = obj_func;

// initialize settings
solver_settings_type settings = solver_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_macro = 100;
strcpy(settings.subsystem_solver, "tcg");

// hand the callback functions whatever the host needs to reach its own data
settings.context = &host_data;

// run solver
c_int error = solver(update_orbs_funptr, obj_func_funptr, n_param, &settings);

// read back output fields
printf("Number of orbital updates: %lld\n", (long long)settings.n_update_orbs);
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol` and `n_macro`) can then be overridden.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Finally, the `solver` is called with a pointer to the initialized settings and callback functions and directly returns an error code in typical C fashion.
- After the call, output fields on `settings` (here, `n_update_orbs`) are populated and can be read like any other attribute.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, solver

# initialize settings
settings = SolverSettings()

# override default settings
settings.conv_tol = 1e-6
settings.n_macro = 100
settings.subsystem_solver = "tcg"

# run solver
solver(update_orbs, obj_func, n_param, settings)

# read back output fields
print(f"Number of orbital updates: {settings.n_update_orbs}")
```

- Callback functions (`update_orbs`, `obj_func`) are defined elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, and `n_macro`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- Finally, the `solver` is called with the initialized settings and callback functions and errors can be caught in pythonic fashion in the form of a `RuntimeError`.
- After the call, output fields on `settings` (here, `n_update_orbs`) are populated and can be read like any other attribute.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for optimization using non-redundant parameters. When this is used, all other passed routines (`update_orbs`, `hess_x`, and `precond`) must be self-projecting. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`subsystem_solver`** (string): Specifies which subsystem solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"tcg"`: truncated conjugate gradient method.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations.
- **`start_trust_radius`** (real): Initial trust radius.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region.
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region.
- **`verbose`** (integer): Controls the verbosity of output during optimization. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it and never keeps it past the call, so it only has to stay valid for the duration of the call. Two solves can therefore run at the same time, or be nested inside one another, as long as each is given its own settings object. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).
- **`stability_settings`** (stability_settings_type): Settings object controlling the internal stability check that is automatically performed upon convergence when `stability` is `True` or when starting at a stationary point (see the Stability Check section below). If `stability_settings%precond`, `stability_settings%project`, `stability_settings%logger`, or `stability_settings%context` are left unset, they default to the corresponding `precond`, `project`, `logger`, and `context` supplied to `solver`. The internal stability check hands its own context to every callback function it calls, so when `stability_settings%context` is set, the Hessian linear transformation returned by `update_orbs` and any inherited `precond`, `project` or `logger` receive it instead of the solver's `context` and must accept it. Leaving it unset keeps the solver's `context` everywhere. `stability_settings%verbose` is raised to at least the solver's own `verbose` level.

### Output
After `solver` returns, the following fields on the settings object have been populated and can be read by the caller:

- **`n_update_orbs`** (integer): Total number of orbital update calls (`update_orbs`) performed by this `solver` call.
- **`n_hess_x`** (integer): Total number of Hessian linear transformations (`hess_x`) performed by this `solver` call.
- **`stability_settings%n_hess_x`** (integer): Number of Hessian linear transformations performed by the most recent internal stability check of this `solver` call.
- **`max_precision_reached`** (boolean): `true` if the solver stopped because the objective function stopped changing, or the trust radius collapsed, at the limit of floating-point precision, rather than because the RMS gradient dropped below `conv_tol` or a custom `conv_check` reported convergence. `error` is still `0` and the returned point is still the best one found; this flag lets the host program decide whether to accept it as-is, warn, or retighten and retry with a different starting point or subsystem solver.

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. If not, it returns a boolean indicating instability and optionally, writes the eigenvector corresponding to the negative eigenvalue in-place to the provided memory.

### Required Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`stable`** (boolean): Returns whether the current point is stable.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`kappa`** (real array): If the memory is provided, the descent direction along the unstable mode is written in-place in this array when the current point is not stable (as can be checked from `stable`), and zeros are written when it is stable.
- **`settings`** (settings_type): Settings object which controls optional arguments as described below.

---

The following Fortran snippet demonstrates how to use the stability check interface:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, hess_x_type, stability_check

real(rp), allocatable :: h_diag(:), kappa(:)
procedure(hess_x_type), pointer :: hess_x_funptr
integer(ip) :: n_param, error
logical :: stable
type(stability_settings_type) :: settings

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! initialize settings
call settings%init(error)

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_iter = 100
settings%diag_solver = "jacobi-davidson"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run stability check
call stability_check(h_diag, hess_x_funptr, stable, error, settings, kappa=kappa)

! read back output fields
print *, "Number of Hessian linear transformations:", settings%n_hess_x
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Stability settings are initialized via the `init()` method of the derived type and can be overridden (here, `conv_tol` and `n_iter`).
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- The `stable` logical output receives the result of the stability check.
- The descent direction `kappa` is optional and is only returned if provided.
- After the call, output fields on `settings` (here, `n_hess_x`) are populated and can be read like any other component.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdio.h>
#include <string.h>
#include "opentrustregion.h"

c_int n_param;
c_bool stable;

// set callback function pointer to existing implementation
hess_x_fp hess_x_funptr = hess_x;

// initialize settings
stability_settings_type settings = stability_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_iter = 100;
strcpy(settings.diag_solver, "jacobi-davidson");

// hand the callback functions whatever the host needs to reach its own data
settings.context = &host_data;

// pointers to Hessian diagonal and descent direction
c_real* h_diag;
c_real* kappa;

// run stability check
c_int error = stability_check(h_diag, hess_x_funptr, n_param, &stable, &settings, kappa);

// read back output fields
printf("Number of Hessian linear transformations: %lld\n", (long long)settings.n_hess_x);
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol` and `n_iter`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- The descent direction `kappa` can be defined elsewhere if needed; otherwise, it can be set to `NULL`.
- After the call, output fields on `settings` (here, `n_hess_x`) are populated and can be read like any other component.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import StabilitySettings, stability_check

# initialize settings
settings = StabilitySettings()

# override default settings
settings.conv_tol = 1e-6
settings.n_iter = 100
settings.diag_solver = "jacobi-davidson"

# Hessian diagonal and descent direction arrays
h_diag = np.asarray(h_diag, dtype=np.float64)
kappa = np.empty(n_param, dtype=np.float64)

# run stability check
stable = stability_check(h_diag, hess_x, n_param, settings, kappa=kappa)

# read back output fields
print(f"Number of Hessian linear transformations: {settings.n_hess_x}")
```

- `hess_x` is an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeError`.
- The descent direction `kappa` is optional and is only returned if provided.
- After the call, output fields on `settings` (here, `n_hess_x`) are populated and can be read like any other attribute.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for stability check using non-redundant parameters. When this is used, all other passed routines (`hess_x` and `precond`) must be self-projecting. Receives the host context as its last argument.
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method that switches to Jacobi-Davidson after `jacobi_davidson_start` iterations.
- **`conv_tol`** (real): Convergence criterion for the residual norm.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of iterations after which the diagonalization solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it and never keeps it past the call, so it only has to stay valid for the duration of the call. Two stability checks can therefore run at the same time, or be nested inside one another, as long as each is given its own settings object. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).

### Output
After `stability_check` returns, the following field on the settings object has been populated and can be read by the caller:

- **`n_hess_x`** (integer): Total number of Hessian linear transformations (`hess_x`) performed by this `stability_check` call.

## Error Code Structure

The library uses structured integer return codes to indicate whether a function has encountered an error. These codes follow the format **`OOEE`**, where:

- **`OO`** = Origin of the error (which component/function reported the error)
- **`EE`** = Specific error code

### General Rules

- A return code of `0` means success.
- Return codes between `1` and `99` are currently unused.
- All current error codes start from `100` and follow the `OOEE` structure.

### Origins (`OO`)

| Code Prefix (`OO`) | Component           |
|--------------------|---------------------|
| `01`               | `solver`            |
| `02`               | `stability_check`   |
| `11`               | `obj_func`          |
| `12`               | `update_orbs`       |
| `13`               | `hess_x`            |
| `14`               | `precond`           |
| `15`               | `conv_check`        |
| `16`               | `project`           |

### Error Codes (`EE`)

The error field (`EE`) is `01` for a general, unspecified error. Some origins define additional, more specific codes where the host program can plausibly react differently to them (e.g. by adjusting a setting and retrying):

| Error Code | Meaning                                                                |
|------------|------------------------------------------------------------------------|
| `0101`     | General error in `solver` |
| `0102`     | Orbital optimization did not converge within the maximum number of macro iterations (`n_macro`) |
| `0201`     | General error in `stability_check` |
| `0202`     | Stability check did not converge within the maximum number of iterations (`n_iter`) |

Future versions may define more specific codes for other actionable failure modes.

### Example Error Codes

| Error Code | Meaning                   |
|------------|---------------------------|
| `0101`     | General error in `solver` |
| `1201`     | Error in `update_orbs`    |

## AI Usage Disclosure

Recent development of OpenTrustRegion has been assisted by AI coding agents. They have been used for code generation, refactoring, writing tests and drafting documentation.

All AI-assisted contributions are reviewed, edited where necessary and validated by the developers before they are merged. The developers make all design decisions, in particular those concerning the underlying methodology, and hold AI-assisted code to the same standards as any other code in the library.
