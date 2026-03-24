# 2D Finite Difference Heat Conduction Solver

This project solves 2D steady-state heat conduction problems using the Finite Difference Method (FDM). It includes different numerical methods, supports common boundary conditions, and provides graphical output of the results.

---

##  Overview

The solver computes the temperature distribution on a rectangular plate by numerically solving the steady heat conduction equation.

### Pure Conduction (Laplace Equation)
∇²T = 0

### Internal Heat Generation (Poisson Equation)
∇²T + q/k = 0

Where:
* T = Temperature
* q = Heat generation per unit volume
* k = Thermal conductivity

The solution is obtained by iteratively updating grid values until the specified convergence tolerance is met.

---

##  Repository Structure

The codebase is organized into a modular pipeline, separating grid generation, boundary logic, numerical solvers, and data export.

.
├── main.py
├── mesh.py
├── boundary.py
├── postprocess.py
└── solver/
    ├── analytical.py
    ├── gauss_seidel.py
    ├── jacobi.py
    ├── sor.py
    └── utils.py

* `main.py`: The core execution script handling user inputs and parallel processing.
* `mesh.py`: Initializes the computational grid and Dirichlet fixed-temperature nodes.
* `boundary.py`: Manages Neumann boundary conditions (convective layers) and calculates boundary states.
* `postprocess.py`: Handles Matplotlib visualization (contours, isotherms) and data export.
* `solver/`
    * `analytical.py`: Computes the exact Fourier series solution for validation.
    * `gauss_seidel.py`: Gauss-Seidel algorithm (includes internal heat generation support).
    * `jacobi.py`: Standard Jacobi iterative scheme.
    * `sor.py`: Successive Over-Relaxation (SOR) with a custom relaxation factor.
    * `utils.py`: Helper functions for error calculations and Biot number extraction.

---

##  Key Features

* **Parallel Solver Racing:** Uses Python's `multiprocessing` to run Jacobi, Gauss-Seidel, and SOR simultaneously, automatically identifying and returning the fastest converged solution.
* **Advanced Boundary Conditions:** Supports standard fixed temperatures (Dirichlet) as well as convective boundary layers (Neumann) utilizing the Biot number (Bi).
* **Analytical Validation:** Automatically compares the numerical solution against the exact 2D analytical Fourier series.
* **Robust Post-Processing:** Automatically generates and saves high-quality `.png` contour maps, isotherm plots, and exports the raw temperature/error fields to `.csv` and `.txt`.

---

##  Numerical Schemes

| Solver | Characteristics | Best Use Case |
| :--- | :--- | :--- |
| **Jacobi** | Uses old temperature values for all spatial updates. | Stable baseline testing, easy to parallelize natively. |
| **Gauss-Seidel** | Uses the most recently computed spatial values. | Faster convergence than Jacobi with lower memory overhead. |
| **SOR** | Applies a relaxation factor (ω) to accelerate the GS update. | High-resolution grids where rapid convergence is required. |

---

##  Usage

To launch the solver pipeline, run the main execution file:

    python main.py

*(You can also use python3 main.py depending on your environment)*

**Interactive Prompts:**
The terminal will guide you through setting up the simulation domain:
* Domain dimensions (Length and Width)
* Grid resolution (Number of nodes in X and Y)
* Wall temperatures (Top, Bottom, Left, Right)
* Convergence tolerance (e.g., 1e-5)
* Convection parameters (Heat transfer coefficient, ambient temp, thermal conductivity)
* Solver selection (Choose a specific solver or race them all)

All visualizations and data exports will be automatically generated and saved in a local `results/` folder.

---

##  Requirements and Dependencies

* `numpy`
* `matplotlib`

Install the required Python libraries using pip:

    pip install numpy matplotlib

---

## Future Work

* Improve convergence tracking.
* Extend the solver to transient (unsteady) conduction.
* Add support for more complex boundary conditions.

---

## Author
**Ruthvik**  

Mechanical Engineering Student