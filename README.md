![License](https://img.shields.io/github/license/eriksen-lab/opentrustregion)
![CI](https://github.com/eriksen-lab/opentrustregion/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/github/eriksen-lab/opentrustregion/graph/badge.svg?token=NJIM6FDADD)](https://codecov.io/github/eriksen-lab/opentrustregion)
[![DOI](https://zenodo.org/badge/855894860.svg)](https://doi.org/10.5281/zenodo.17142554)

# OpenTrustRegion: A Reusable Library for Second-Order Trust Region Orbital Optimization

This library provides a robust and flexible implementation for second-order trust region orbital optimization, with extensive customization options to suit various use cases.

The following paper documents the theory and implementation of the methodology in OpenTrustRegion, and should be cited in any work using OpenTrustRegion:  

- A Reusable Library for Second-Order Orbital Optimization Using the Trust Region Method\
  Greiner, J.; Høyvik, I.-M.; Lehtola, S.; Eriksen, J. J.\
  [arXiv:2509.13931](https://arxiv.org/abs/2509.13931)

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
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
- **`obj_func`** (function):  
  Accepts and applies a variable update (e.g., orbital rotation) and returns:
  - Objective function value (real)
  - An integer error code (0 for success, positive integers < 100 for errors)
- **`n_param`** (integer): Specifies the number of parameters to be optimized.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.

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

! run solver
call solver(update_orbs_funptr, obj_func_funptr, n_param, error, settings)
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized using the `init()` method of the derived type, and default settings can be overridden (here, `conv_tol` and `n_macro`).
- Finally, the `solver` is called with the initialized settings and callback functions.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdint.h>
#include <string.h>
#include "opentrustregion.h"

c_ip n_param;

// set callback function pointers to existing implementations
update_orbs_type update_orbs_funptr = (void*)update_orbs;
obj_func_type obj_func_funptr = (void*)obj_func;

// initialize settings
struct solver_settings_type settings
struct solver_settings_type settings = solver_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_macro = 100;
strcpy(settings.subsystem_solver, "tcg");

// run solver
c_ip error = solver(update_orbs_funptr, obj_func_funptr, n_param, settings);
```

- Callback function pointers (`update_orbs_funptr`, `obj_func_funptr`) point to existing implementations elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via a small helper function `solver_settings_init()`, which returns a struct with default values. Individual settings (here, `conv_tol` and `n_macro`) can then be overridden.
- Finally, the `solver` is called with the initialized settings and callback functions and directly returns an error code in typical C fashion.

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
```

- Callback functions (`update_orbs`, `obj_func`) are defined elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Solver settings are initialized via the `SolverSettings` class, which returns an object with default values; individual settings (here, `conv_tol`, and `n_macro`) can then be overridden.
- Finally, the `solver` is called with the initialized settings and callback functions and errors can be caught in pythonic fashion in the form of a `RuntimeException`.

### Optional Settings
The optimization process can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors).
- **`conv_check`** (function): Returns whether the optimization has converged due to some supplied convergence criterion. Additionally, outputs an integer code indicating the success or failure of the function, positive integers less than 100 represent error conditions.
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
- **`verbose`** (integer): Controls the verbosity of output during optimization.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (subroutine): Accepts a log message. Logging is otherwise routed to stdout.

## Stability Check
A separate `stability_check` subroutine is available to verify whether the current solution corresponds to a minimum. If not, it returns a boolean indicating instability and optionally, writes the eigenvector corresponding to the negative eigenvalue in-place to the provided memory.

### Required Arguments

- **`h_diag`** (real array): Represents the Hessian diagonal at the current point.
- **`hess_x`** (subroutine): Performs Hessian-vector products at the current point:
  - Accepts a trial vector and writes the result of the Hessian transformation into an output array (real array, written in-place)
  - Returns an integer error code (0 for success, positive integers < 100 for errors)
- **`stable`** (boolean): Returns whether the current point is stable.
- **`error`** (integer): An integer code indicating the success or failure of the solver. The error code structure is explained below.
- **`kappa`** (real array): If the memory is provided and the current point is not stable (as can be checked from return code of `stable`), the descent direction is written in-place in this array.

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

! run stability check
call stability_check(h_diag, hess_x_funptr, n_param, stable, error, settings, kappa=kappa)
```

- `hess_x_funptr` points to an existing Hessian–vector product implementation elsewhere in the program.
- `n_param` is also assumed to be defined elsewhere.
- Stability settings are initialized via the `init()` method of the derived type and can be overridden (here, `conv_tol` and `n_iter`).
- The `stable` logical output receives the result of the stability check.
- The descent direction `kappa` is optional and is only returned if provided.

---

The following C snippet demonstrates the equivalent usage through the C interface:

```c
#include <stdint.h>
#include <string.h>
#include "opentrustregion.h"

c_ip n_param;
c_bool stable;

// set callback function pointer to existing implementation
hess_x_type hess_x_funptr = hess_x;

// initialize settings
struct stability_settings_type settings = stability_settings_init();

// override default settings
settings.conv_tol = 1e-6;
settings.n_iter = 100;
strcpy(settings.diag_solver, "jacobi-davidson");

// pointers to Hessian diagonal and descent direction
double* h_diag;
double* kappa;

// run stability check
c_ip = stability_check(h_diag, hess_x_funptr, n_param, &stable, settings, kappa);
```

- `hess_x_funptr` points to an existing Hessian–vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via a small helper function `stability_settings_init()`, which returns a struct with default values; individual settings (here, `conv_tol` and `n_iter`) can then be overridden.
- The `stable` output receives the result of the stability check which directly returns an error code in typical C fashion.
- The descent direction `kappa` can be defined elsewhere if needed; otherwise, it can be set to `nullptr`.

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
```

- `hess_x` is an existing Hessian–vector product implementation elsewhere in the program.
- `n_param` and `h_diag` are assumed to be defined elsewhere.
- Stability settings are initialized via the `StabilitySettings` class, which returns an object with default values; individual settings (here, `conv_tol`) can then be overridden.
- The `stable` output receives the result of the stability check and errors can be caught in pythonic fashion in the form of a `RuntimeException`.
- The descent direction `kappa` is optional and is only returned if provided.

### Optional Settings
The stability check can be fine-tuned using the following settings:

- **`precond`** (subroutine): Applies a preconditioner to a residual vector. Writes the result in-place into a provided array and returns an integer error code (0 for success, positive integers < 100 for errors).
- **`diag_solver`** (string): Specifies which diagonalization solver to use. Options include:
  - `"davidson"`: standard Davidson method,
  - `"jacobi-davidson"`: Davidson method with fallback to Jacobi-Davidson if convergence is difficult, or automatically after `jacobi_davidson_start` micro iterations.
- **`conv_tol`** (real): Convergence criterion for the residual norm.
- **`n_random_trial_vectors`** (integer): Number of random trial vectors used to start the Davidson iterations.
- **`n_iter`** (integer): Maximum number of Davidson iterations.
- **`jacobi_davidson_start`** (integer): Number of micro iterations after which the subsystem solver switches to the Jacobi-Davidson method.
- **`verbose`** (integer): Controls the verbosity of output during the stability check.
- **`seed`** (integer): Seed value for generating random trial vectors.
- **`logger`** (function): Accepts a log message. Logging is otherwise routed to stdout.

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

### Error Codes (`EE`)

The error field (`EE`) is currently always set to `01`. Future versions may define more specific codes for different failure modes.

### Example Error Codes

| Error Code | Meaning                |
|------------|------------------------|
| `0101`     | Error in `solver`      |
| `1201`     | Error in `update_orbs` |

## Program Interfaces

The PySCF interface is available as an extension hosted at https://github.com/eriksen-lab/pyscf_opentrustregion. To install it, simply add its path to the **`PYSCF_EXT_PATH`** environment variable:
```sh
export PYSCF_EXT_PATH=path/to/pyscf_opentrustregion
```
Usage examples can be found in the **`examples`** directory of the PySCF interface repository.

The interface supports Hartree–Fock and DFT calculations via the **`mf_to_otr`** function, which wraps PySCF **`HF`** and **`KS`** objects into their OpenTrustRegion counterparts. Similarly, localization methods are available through the **`BoysOTR`**, **`PipekMezeyOTR`**, and **`EdmistonRuedenbergOTR`** classes, and state-specific CASSCF calculations are supported via the **`casscf_to_otr`** function applied to a PySCF **`CASSCF`** object. All returned objects are fully compatible with the original PySCF classes and can be used interchangeably.

Optional settings can be adjusted by modifying object attributes directly. Orbital optimization and internal stability analysis are performed using the **`kernel`** and **`stability_check`** member functions, respectively.
