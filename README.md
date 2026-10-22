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

The library needs a Fortran 2008 compiler, a C compiler, CMake 3.24 or newer and BLAS and LAPACK. It is tested with gfortran 10 and newer, Intel ifx and LLVM Flang 21 and newer, and against OpenBLAS, the reference BLAS and LAPACK, Intel MKL and Apple Accelerate.

### Fortran or C Usage
To build the library for Fortran or C:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The installation can be tested by running the ```testsuite.py``` file in the ```pyopentrustregion``` directory.

A Fortran or C program has to be rebuilt against the modules or the header of the library version it links, since the settings and state types can change from one version to the next.

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
| **OpenTrustRegion_BUILD_TESTING** | `BOOL` | `ON` when built as the top-level project, `OFF` otherwise | Build the project’s testsuite. |
| **OpenTrustRegion_INSTALL_TESTING_LIBRARY** | `BOOL` | `OFF` | Install the library the testsuite runs (`libotrtestsuite`) together with the library. |
| **OpenTrustRegion_INSTALL_CMAKEDIR** | `STRING` | (auto) | Project install directory. |
| **CMAKE_BUILD_TYPE** | `STRING` | `Release` | Choose the build type (`Debug`, `Release`, etc.). |
| **INTEGER_SIZE** | `STRING` | *(auto)* | Set the integer precision to `4` (32-bit) or `8` (64-bit). Required when providing custom BLAS/LAPACK libraries. Otherwise defaults to 32-bit integers and tries to locate compatible BLAS and LAPACK libraries. Falls back to 64-bit integers if 32-bit libraries cannot be found. The resulting library name reflects the chosen integer precision (`libopentrustregion_32.*` or `libopentrustregion_64.*`). C programs compiled against the header of a library with 64-bit integers need `OTR_ILP64` defined, which linking the CMake target `OpenTrustRegion::opentrustregion` does automatically. |
| **BLAS_LIBRARIES** | `PATH` | *(auto)* | Path(s) to BLAS libraries. If not provided, CMake attempts to locate a suitable BLAS automatically. |
| **LAPACK_LIBRARIES** | `PATH` | *(auto)* | Path(s) to LAPACK libraries. If not provided, CMake attempts to locate a suitable LAPACK automatically. |
| **OpenTrustRegion_HOST_PROVIDES_BLAS** | `BOOL` | `OFF` | When enabled, OpenTrustRegion will not attempt to detect or link BLAS/LAPACK and the testsuite is automatically disabled. The calling program must provide BLAS/LAPACK routines that expose the unsuffixed symbol names (for example `ddot`, `dsyev`) with an integer width matching `INTEGER_SIZE`; otherwise linking will fail loudly. |
| **OpenTrustRegion_ENABLE_XHOST** | `BOOL` | `ON` | Optimize the Release build for the current machine's instruction set (`-march=native` for GNU, `-xHost` for Intel). Automatically disabled when cross-compiling, regardless of this setting. Turn off when building for a different machine than the one compiling (e.g. packaging/conda builds). |

### Nested and Concurrent Calls

The library keeps no state of its own between calls, so a solver or stability check call can be nested inside a callback function of another call, or run concurrently with another call in a different thread, as long as each call is given its own state record and, in Python, its own settings object, into which every call writes the interfaces to its callback functions. A nested call re-enters procedures of the library that are still running, which Fortran only permits for procedures compiled for recursion, and a concurrent call needs their local variables on its own thread's stack. The library then has to be built with `-frecursive` for gfortran (implied by `-fopenmp`) or `-recursive` for Intel:

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

The optimization process is initiated by calling a `solver` subroutine. This routine takes the following arguments:

### Arguments

- **`update_orbs`** (subroutine):  
  Accepts and applies a variable update (e.g., orbital rotation), updates the host's internal state, and provides:
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
- **`settings`** (settings_type, optional): Settings object which controls optional arguments as described below. The solver only reads it. Without it, the default settings are used.
- **`state`** (solver_state_type, optional): Record that the solver keeps up to date during the optimization and that holds the final state of the optimization once the solver returns, as described below.

---

The following Fortran snippet demonstrates how to use the `solver` interface:

```fortran
use opentrustregion, only: ip, rp, update_orbs_type, obj_func_type, &
                           solver_settings_type, solver_state_type, solver

procedure(update_orbs_type), pointer :: update_orbs_funptr
procedure(obj_func_type), pointer :: obj_func_funptr
integer(ip) :: n_param, error
type(solver_settings_type) :: settings
type(solver_state_type) :: state

! set callback function pointers to existing implementations
update_orbs_funptr => update_orbs
obj_func_funptr => obj_func

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_macro = 100
settings%subsystem_solver = "tcg"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run solver
call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings, state)

! read the final state
print *, "Number of orbital updates:", state%n_update_orbs
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- A declared settings object holds the default settings, which can be overridden (here, `conv_tol` and `n_macro`); assigning `solver_settings_type()` resets an object to them.
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- Finally, the `solver` is called with the callback functions, the settings and a state record. Both of these are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the optimization (here, the number of orbital updates `n_update_orbs`), see [Output](#output). The settings are left unchanged.

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
solver_state_type state;

// override default settings
settings.conv_tol = 1e-6;
settings.n_macro = 100;
strcpy(settings.subsystem_solver, "tcg");

// hand the callback functions whatever the host needs to reach its own data
settings.context = &host_data;

// run solver
c_int error = solver(update_orbs_funptr, obj_func_funptr, n_param, &settings, &state);

// read the final state
printf("Number of orbital updates: %lld\n", (long long)state.n_update_orbs);
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol` and `n_macro`) can then be overridden.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Finally, the `solver` is called with the callback functions and pointers to the settings and to a state struct, and directly returns an error code in typical C fashion. Either pointer may be `NULL`: without settings, the default settings are used, and without a state struct, none is returned. Since the host context is a setting, a host that needs one passes settings.
- After the call, `state` holds the final state of the optimization (here, `n_update_orbs`), also after an error, see [Output](#output). The state struct needs no initialization, and the settings are left unchanged.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import SolverSettings, SolverState, solver

# initialize settings and state
settings = SolverSettings()
state = SolverState()

# override default settings
settings.conv_tol = 1e-6
settings.n_macro = 100
settings.subsystem_solver = "tcg"

# run solver
solver(update_orbs, obj_func, n_param, settings, state)

# read the final state
print(f"Number of orbital updates: {state.n_update_orbs}")
```

- Callback functions (`update_orbs`, `obj_func`) are defined elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, and `n_macro`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- Finally, the `solver` is called with the callback functions, the settings and a state object, and errors can be caught in pythonic fashion in the form of a `RuntimeError`. The settings and the state are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the optimization (here, `n_update_orbs`), see [Output](#output), whose fields are read-only attributes. The settings are left unchanged.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for optimization using non-redundant parameters. When this is used, all other passed routines (`update_orbs`, `hess_x`, and `precond`) must be self-projecting. Receives the host context as its last argument.
- **`get_extra_trial_vectors`** (subroutine): Returns extra trial vectors that seed the leading part of the initial trial space, written in-place to the provided matrix, whose column count is `n_extra_trial_vectors` and is passed to the callback itself by the C and Python interfaces, since those only receive a flat pointer. In Python, the matrix holds one trial vector per row instead. The vectors do not need to be orthonormalized or projected, and a vanishing column is dropped, which is how the callback signals that it has no vector to contribute for that slot. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion, which replaces the criterion on the RMS gradient (`conv_tol`). A point the solver reaches by a step along an unstable mode does not count as converged, whatever the check returns. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`stability`** (boolean): Determines whether a stability check is performed upon convergence. The solver checks the stability of a starting point that is already converged and of a point with an exactly vanishing gradient in any case. From an unstable point, it steps along the unstable mode until the gradient RMS of the quadratic model is ten times the larger of the current one and `conv_tol`, or by the trust radius if both vanish or the step would be longer, and does not take the point it reaches as converged.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians. A non-symmetric Hessian is only supported by the `"davidson_ls"` subsystem solver. The internal stability check, which diagonalizes the same Hessian, inherits a non-symmetric Hessian, i.e. this setting switched off also switches off `stability_settings%hess_symm`.
- **`line_search`** (boolean): Determines whether a line search is performed after every macro iteration.
- **`grad_noise`** (boolean): Determines whether the gradient is perturbed by a random vector before every subsystem solve of the Davidson-family subsystem solvers. The perturbation keeps the solution from remaining exactly within a symmetry subspace of the starting point, which a Hessian with no negative curvature outside that subspace, or no information about it, cannot otherwise leave. Its norm is the micro-iteration reduction factor (`local_red_factor` or `global_red_factor`) of the region the previous subproblem was solved in, the first counting as global, times the gradient norm, so that it goes to zero together with the gradient at convergence; the micro iterations are never asked to reduce the residual below this relative size, so that they do not resolve the perturbation itself, whichever region the subproblem turns out to be solved in.
- **`subsystem_solver`** (string): Specifies which subsystem solver to use. Options include:
  - `"davidson_ls"`: generalized Davidson method applied to the linear system,
  - `"jacobi-davidson_ls"`: generalized Davidson method applied to the linear system with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"davidson_ah"`: Davidson method applied to the diagonalization of the augmented Hessian,
  - `"jacobi-davidson_ah"`: Davidson method applied to the diagonalization of the augmented Hessian with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations,
  - `"tcg"`: truncated conjugate gradient method.
- **`conv_tol`** (real): Specifies the convergence criterion for the RMS gradient, which a custom `conv_check` replaces.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to initialize the micro iterations of the Davidson-family subsystem solvers.
- **`n_extra_trial_vectors`** (integer): Number of non-random trial vectors added alongside the gradient direction to initialize the micro iterations, on top of the `n_random_trial_vectors` random ones. These are taken from `get_extra_trial_vectors` when that callback is provided, and are otherwise unit vectors along the lowest Hessian diagonal elements, added only where those indicate negative curvature. Vectors that vanish or are linearly dependent on the gradient direction are dropped.
- **`start_trust_radius`** (real): Initial trust radius.
- **`n_macro`** (integer): Maximum number of macro iterations.
- **`n_micro`** (integer): Maximum number of micro iterations per attempted trust-region step.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`global_red_factor`** (real): Reduction factor for the residual during micro iterations in the global region, also the relative size of the gradient perturbation there (see `grad_noise`).
- **`local_red_factor`** (real): Reduction factor for the residual during micro iterations in the local region, also the relative size of the gradient perturbation there (see `grad_noise`).
- **`verbose`** (integer): Controls the verbosity of output during optimization. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors and the random gradient perturbation.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it and never keeps it past the call, so it only has to stay valid for the duration of the call. Two solves can therefore run at the same time, or be nested inside one another, as long as each is given its own state record. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).
- **`stability_settings`** (stability_settings_type): Settings object controlling the internal stability check that is automatically performed upon convergence when `stability` is `True` and in any case at a starting point that is already converged and at a point with an exactly vanishing gradient (see the Stability Check section below). If `stability_settings%precond`, `stability_settings%project`, `stability_settings%logger`, or `stability_settings%context` are left unset, they default to the corresponding `precond`, `project`, `logger`, and `context` supplied to `solver`. The internal stability check hands its own context to every callback function it calls, so when `stability_settings%context` is set, the Hessian linear transformation returned by `update_orbs` and any inherited `precond`, `project` or `logger` receive it instead of the solver's `context` and must accept it. Leaving it unset keeps the solver's `context` everywhere. `stability_settings%verbose` is raised to at least the solver's own `verbose` level, and `stability_settings%hess_symm` is switched off together with the solver's `hess_symm`. The solver checks the nested settings, after these defaults are resolved, together with its own before the optimization starts, also without `stability`, since the internal stability check runs in any case at the points named above.

### Output
The solver keeps the state it is given up to date during the optimization, starting from the default values on every call, and leaves the final state of the optimization in it when it returns, also after an error. A callback function can therefore read the current state of the optimization: in Fortran, the state record is declared with the `target` attribute and reached through a pointer, for example one in the host context; in C, the state is copied into the host's struct before every callback function; and in Python, the state object can be read from within a callback function. The state holds the following fields:

- **`macro_iteration`** (integer): Macro iteration of the current point, 0 for the starting point.
- **`func`** (real): Objective function value at the current point.
- **`grad_rms`** (real): RMS gradient at the current point.
- **`trust_radius`** (real): Current trust radius.
- **`micro_iteration`** (integer): Number of micro iterations spent on the trust-region subproblem being solved or, until the next one starts, on the last one, including those of rejected steps, and 0 after a step along an unstable mode.
- **`jacobi_davidson_iteration`** (integer): Number of these micro iterations that were Jacobi-Davidson iterations.
- **`relative_residual`** (real): Euclidean norm of the residual of the trust-region subproblem being solved, relative to that of the gradient it is solved for, which `grad_noise` perturbs, for the Davidson-family subsystem solvers and to that of the initial residual, the gradient with a small random perturbation, for `"tcg"`, whose residual is that of its last iterate inside the trust region, and 0 after a step along an unstable mode.
- **`level_shift`** (real): Level shift of the last step, which the Davidson-family subsystem solvers apply to the Hessian to keep a step within the trust region, 0 for a Newton step inside the trust region, for a step along an unstable mode and for `"tcg"`.
- **`step_size`** (real): Size of the last step in the norm of the trust region, after the line search if one is performed. A step along an unstable mode is measured in the Euclidean norm, also for `"tcg"`.
- **`n_update_orbs`** (integer): Number of orbital updates (`update_orbs`).
- **`n_obj_func`** (integer): Number of objective function evaluations (`obj_func`), including those of the line search.
- **`n_hess_x`** (integer): Number of Hessian linear transformations (`hess_x`), those of an internal stability check added once it returns.
- **`max_precision_reached`** (boolean): `true` once neither the actual nor the predicted change of the objective function in a trust-region step can be resolved at floating-point precision relative to the objective function value, once a step is rejected whose predicted change can not be resolved, as when the objective function and the orbital update disagree by more than that, or once the gradient vanishes exactly without passing the convergence check, which leaves the subsystem solvers no subproblem to solve. This ends the optimization like convergence, without an error of its own, and the returned point is the best one found; this flag lets the host program decide whether to accept it as-is, warn, or retighten and retry with a different starting point or subsystem solver. It turns `false` again if the internal stability check, which runs with `stability` set and at an exactly vanishing gradient, then finds the point unstable and the solver continues along the unstable mode, or if the convergence check (`conv_check`, or without one the criterion on the RMS gradient) passes, which takes precedence.
- **`stability_state`** (stability_state_type): State of the most recent internal stability check, see the [stability check's output](#output-1).

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. It returns whether the point is stable and, optionally, the lowest Hessian eigenvalue and the corresponding eigenvector, which at an unstable point is the direction along which the objective function decreases. It takes the following arguments:

### Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
  - Receives the host context as its last argument.
- **`stable`** (boolean): Returns whether the current point is stable, i.e. whether the lowest Hessian eigenvalue is above -1e-2; a less negative eigenvalue is not distinguished from numerical noise, such as that of the integration grid of a density functional.
- **`error`** (integer): An integer code indicating the success or failure of the stability check. The error code structure is explained below.
- **`settings`** (settings_type, optional): Settings object which controls optional arguments as described below. The stability check only reads it. Without it, the default settings are used.
- **`state`** (stability_state_type, optional): Record that the stability check keeps up to date during its iterations and that holds its final state once it returns, as described below.
- **`kappa`** (real array, optional): If provided, receives the normalized eigenvector of the lowest Hessian eigenvalue, written in-place, whether the current point is stable or not. At an unstable point (as can be checked from `stable`), this is the descent direction along the unstable mode.

---

The following Fortran snippet demonstrates how to use the stability check interface:

```fortran
use opentrustregion, only: ip, rp, stability_settings_type, stability_state_type, &
                           hess_x_type, stability_check

real(rp), allocatable :: h_diag(:), kappa(:)
procedure(hess_x_type), pointer :: hess_x_funptr
integer(ip) :: n_param, error
logical :: stable
type(stability_settings_type) :: settings
type(stability_state_type) :: state

! set callback function pointer to existing implementation
hess_x_funptr => hess_x

! override default settings
settings%conv_tol = 1e-6_rp
settings%n_iter = 100
settings%diag_solver = "jacobi-davidson"

! hand the callback functions whatever the host needs to reach its own data
settings%context => host_data

! run stability check
call stability_check(h_diag, hess_x_funptr, stable, error, settings, state, &
                     kappa=kappa)

! read the final state
print *, "Lowest Hessian eigenvalue:", state%min_eigval
print *, "Number of Hessian linear transformations:", state%n_hess_x
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- A declared settings object holds the default settings, which can be overridden (here, `conv_tol` and `n_iter`); assigning `stability_settings_type()` resets an object to them.
- `host_data` is a variable of any type declared with the `target` attribute. Every callback function receives it as its last argument, declared as `class(*), intent(in), pointer :: context`, and recovers its own type with `select type`. The `intent(in)` applies to the pointer rather than the data, so a callback function may modify `host_data` but not point `context` elsewhere. Leaving `settings%context` unset is fine; the callback functions then receive an unassociated pointer.
- The `stable` logical output receives the result of the stability check.
- The settings, the state and the lowest eigenvector `kappa` are optional: without settings, the default settings are used, and the others are only returned if provided.
- After the call, `state` holds the final state of the stability check (here, the lowest eigenvalue `min_eigval` and `n_hess_x`), see [Output](#output-1). The settings are left unchanged.

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
stability_state_type state;

// override default settings
settings.conv_tol = 1e-6;
settings.n_iter = 100;
strcpy(settings.diag_solver, "jacobi-davidson");

// hand the callback functions whatever the host needs to reach its own data
settings.context = &host_data;

// pointers to Hessian diagonal and lowest eigenvector
c_real* h_diag;
c_real* kappa;

// run stability check
c_int error = stability_check(h_diag, hess_x_funptr, n_param, &stable, &settings,
                              &state, kappa);

// read the final state
printf("Lowest Hessian eigenvalue: %f\n", state.min_eigval);
printf("Number of Hessian linear transformations: %lld\n", (long long)state.n_hess_x);
```

- `hess_x_funptr` points to an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- `host_data` is any host object. Every callback function receives `&host_data` as its `void *context` argument and casts it back. Leaving `settings.context` as `NULL` is fine; the callback functions then receive `NULL`.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol` and `n_iter`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- `kappa` receives the lowest eigenvector and can be `NULL` to skip it. So can `state`, and `settings` can be `NULL` to use the default settings.
- After the call, `state` holds the final state of the stability check (here, the lowest eigenvalue `min_eigval` and `n_hess_x`), also after an error, see [Output](#output-1). The state struct needs no initialization, and the settings are left unchanged.

---

The following Python snippet demonstrates the equivalent usage through the Python interface:

```python
from pyopentrustregion import StabilitySettings, StabilityState, stability_check

# initialize settings and state
settings = StabilitySettings()
state = StabilityState()

# override default settings
settings.conv_tol = 1e-6
settings.n_iter = 100
settings.diag_solver = "jacobi-davidson"

# Hessian diagonal and lowest eigenvector arrays
h_diag = np.asarray(h_diag, dtype=np.float64)
kappa = np.empty(n_param, dtype=np.float64)

# run stability check
stable = stability_check(h_diag, hess_x, n_param, settings, state, kappa=kappa)

# read the final state
print(f"Lowest Hessian eigenvalue: {state.min_eigval}")
print(f"Number of Hessian linear transformations: {state.n_hess_x}")
```

- `hess_x` is an existing Hessian-vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`) can then be overridden.
- There is no `context` setting in Python. Any callable works as a callback function, so a closure, a bound method or a `functools.partial` already carries whatever host data the callback function needs.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeError`.
- The lowest eigenvector is written into `kappa` if it is provided.
- The settings and the state are optional: without settings, the default settings are used, and the state is only needed to read it.
- After the call, `state` holds the final state of the stability check (here, the lowest eigenvalue `min_eigval` and `n_hess_x`), see [Output](#output-1), whose fields are read-only attributes. The settings are left unchanged.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors). Receives the host context as its last argument.
- **`project`** (subroutine): Applies a projection in-place to a provided vector and returns an integer error code (0 for success, positive integers < 100 for errors). Required for stability check using non-redundant parameters. When this is used, all other passed routines (`hess_x` and `precond`) must be self-projecting. Receives the host context as its last argument.
- **`conv_check`** (function): Returns whether the stability check has converged due to some supplied convergence criterion, which replaces the criterion on the RMS residual (`conv_tol`) and `stop_on_instability`. It has the same arguments as the convergence check of the solver and can read the `state` passed to `stability_check`, or, for the internal stability check, from `stability_state` of the state passed to `solver`, so a host that supplies it has to pass a state. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`get_extra_trial_vectors`** (subroutine): Returns extra trial vectors that seed the leading part of the initial trial space, written in-place to the provided matrix, whose column count is `n_extra_trial_vectors` and is passed to the callback itself by the C and Python interfaces, since those only receive a flat pointer. In Python, the matrix holds one trial vector per row instead. The vectors do not need to be orthonormalized or projected, and a vanishing column is dropped, which is how the callback signals that it has no vector to contribute for that slot. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions. Receives the host context as its last argument.
- **`hess_symm`** (boolean): Determines whether the supplied Hessian is symmetric. This is sometimes not the case for approximate Hessians. A non-symmetric Hessian is supported neither by the `"jacobi-davidson"` diagonalization solver nor by `stop_on_instability`.
- **`stop_on_instability`** (boolean): When no `conv_check` is supplied, lets the default convergence policy accept a Ritz pair as soon as its eigenvalue drops below the instability threshold, without waiting for the residual to converge, since for a symmetric Hessian a Ritz value is an upper bound on the true lowest eigenvalue and the sign is then already certain. Not supported for a non-symmetric Hessian (`hess_symm`). Off by default.
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method that switches to Jacobi-Davidson after `jacobi_davidson_start` iterations.
- **`conv_tol`** (real): Convergence criterion for the RMS residual.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_extra_trial_vectors`** (integer): Number of non-random trial vectors used to start the Davidson iterations, on top of the `n_random_trial_vectors` random ones. These non-random vectors are taken from `get_extra_trial_vectors` when that callback is provided, and otherwise are unit vectors along the lowest Hessian diagonal elements. Vectors that vanish or are linearly dependent are dropped, so the trial space holds the remaining ones plus the random ones, at most as many as there are parameters, and the stability check fails if nothing remains.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of iterations after which the diagonalization solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check. Level 0 prints nothing, 1 prints errors, 2 also warnings, 3 also progress information and 4 also debugging information.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout, and error messages to stderr. Receives the host context as its last argument.
- **`context`** (unlimited polymorphic pointer in Fortran, `void *` in C, absent in Python): Opaque host data, handed back unchanged as the last argument of every callback function so that the host does not have to reach its own state through module-level variables. The library never inspects it and never keeps it past the call, so it only has to stay valid for the duration of the call. Two stability checks can therefore run at the same time, or be nested inside one another, as long as each is given its own state record. This requires a suitable build, see [Nested and Concurrent Calls](#nested-and-concurrent-calls).

### Output
The stability check keeps the state it is given up to date during its iterations, starting from the default values on every call, and leaves its final state in it when it returns, also after an error. As for the solver, a callback function can therefore read the current state. The state holds the following fields:

- **`iteration`** (integer): Current iteration of the diagonalization.
- **`min_eigval`** (real): Current estimate of the lowest Hessian eigenvalue, once the check has converged the lowest eigenvalue.
- **`residual_rms`** (real): RMS of the residual of this estimate.
- **`n_hess_x`** (integer): Number of Hessian linear transformations (`hess_x`).

## Error Code Structure

The library uses structured integer return codes to indicate whether a function has encountered an error. These codes follow the format **`OOEE`**, where:

- **`OO`** = Origin of the error (which component/function reported the error)
- **`EE`** = Specific error code

### General Rules

- A return code of `0` means success.
- Return codes between `1` and `99` are currently unused.
- All current error codes start from `100` and follow the `OOEE` structure.

### Origins (`OO`)

| Code Prefix (`OO`) | Component                                                                  |
|--------------------|----------------------------------------------------------------------------|
| `01`               | `solver`                                                                   |
| `02`               | `stability_check`                                                          |
| `11`               | `obj_func`                                                                 |
| `12`               | `update_orbs`                                                              |
| `13`               | `hess_x`                                                                   |
| `14`               | `precond`                                                                  |
| `15`               | `conv_check`, of the solver and of the stability check                     |
| `16`               | `project`                                                                  |
| `17`               | `get_extra_trial_vectors`, of the solver and of the stability check        |

### Error Codes (`EE`)

The error field (`EE`) is `01` for a general, unspecified error. Some origins define additional, more specific codes where the host program can plausibly react differently to them (e.g. by adjusting a setting and retrying):

| Error Code | Meaning                                                                |
|------------|------------------------------------------------------------------------|
| `0101`     | General error in `solver` |
| `0102`     | Orbital optimization did not converge within the maximum number of macro iterations (`n_macro`) |
| `0201`     | General error in `stability_check` |
| `0202`     | Stability check did not converge, within the maximum number of iterations (`n_iter`) or because its trial space could not be expanded any further |

Future versions may define more specific codes for other actionable failure modes.

### Example Error Codes

| Error Code | Meaning                   |
|------------|---------------------------|
| `0101`     | General error in `solver` |
| `1201`     | Error in `update_orbs`    |

## AI Usage Disclosure

Recent development of OpenTrustRegion has been assisted by AI coding agents. They have been used for code generation, refactoring, writing tests and drafting documentation.

All AI-assisted contributions are reviewed, edited where necessary and validated by the developers before they are merged. The developers make all design decisions, in particular those concerning the underlying methodology, and hold AI-assisted code to the same standards as any other code in the library.
