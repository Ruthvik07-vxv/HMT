# 2D Heat Conduction Solver (Finite Difference Method)

A Python implementation of **2‑D steady state heat conduction solvers**
using the **Finite Difference Method (FDM)**.

This repository contains **two independent solvers** designed to solve
different heat transfer problems on a rectangular plate.\
Both solvers use an **iterative Gauss--Seidel style method** to converge
the temperature field to a specified tolerance.

------------------------------------------------------------------------

# Overview

The code numerically solves the **2‑D steady heat conduction equation**.

### Without heat generation

∇²T = 0

### With internal heat generation

∇²T + q/k = 0

Where:

  Symbol   Meaning
  -------- ---------------------------------
  T        Temperature
  q        Heat generation per unit volume
  k        Thermal conductivity

The solver works on a **rectangular grid** and iteratively updates node
temperatures until convergence.

------------------------------------------------------------------------

# Repository Structure

    .
    ├── solver.py
    ├── solver2.py
    └── README.md

## solver.py

Main FDM solver with the following features:

-   Dirichlet boundary conditions
-   Optional convection boundary conditions
-   Analytical solution comparison (Fourier series)
-   Error calculation
-   CSV and TXT export of results
-   Temperature contour plots
-   Isotherm plots

## solver2.py

Extended solver with additional capabilities:

-   Internal heat generation
-   Arbitrary fixed temperature nodes inside the mesh
-   Separate solvers for:
    -   No heat generation
    -   Heat generation
-   Custom geometry support
-   Convergence monitoring
-   Visualization of results

------------------------------------------------------------------------

# Features

## Numerical Solver

-   2‑D Finite Difference discretization
-   Gauss‑Seidel iterative method
-   Convergence based on maximum error tolerance

## Boundary Conditions

Supported boundary conditions include:

-   Fixed temperature (Dirichlet)
-   Convective boundary conditions
-   Internal fixed temperature nodes

## Heat Generation

The second solver supports volumetric heat generation.

## Analytical Solution Comparison

For pure conduction cases, the solver compares numerical results against
the **analytical Fourier series solution**.

Reported metrics include:

-   Maximum error
-   Average error
-   Node‑by‑node comparison

## Visualization

The solver automatically generates:

-   Temperature contour plots
-   Isotherm lines
-   Combined visualization

Plots are saved as:

-   Temperature Contours.png
-   Temperature Isotherms.png
-   Combined Plot.png

## Data Export

Temperature fields can be exported to:

-   ActualTemperature.txt
-   TemperatureField.csv
-   TheoreticalT.txt
-   TheoreticalT.csv

These files can be analyzed using tools like **Excel, MATLAB, or
Python**.

------------------------------------------------------------------------

# How the Solver Works

1.  A mesh is created with specified boundary temperatures.

2.  Interior nodes are updated using the finite difference formulation:

    T(i,j) = 1/4 \[T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1)\]

3.  Iterations continue until:

    max(\|T_new − T_old\|) \< tolerance

4.  Results are visualized and exported.

------------------------------------------------------------------------

# Requirements

Install required Python libraries:

    pip install numpy matplotlib

------------------------------------------------------------------------

# Running the Solvers

## Run solver.py

    python solver.py

You will be prompted for:

-   Plate dimensions
-   Grid resolution
-   Boundary temperatures
-   Convergence tolerance
-   Optional convection parameters

------------------------------------------------------------------------

## Run solver2.py

    python solver2.py

Choose the case:

    1 → Steady conduction (no heat generation)
    2 → Conduction with heat generation

Then provide:

-   Grid resolution
-   Boundary temperatures
-   Plate dimensions
-   Heat generation rate (if applicable)
-   Thermal conductivity
-   Convergence tolerance

Optional internal fixed temperature nodes can also be specified.

------------------------------------------------------------------------

# Numerical Method

The solver uses:

-   Finite Difference Method (FDM)
-   Gauss‑Seidel iteration
-   Structured rectangular grid

The approach is simple, stable, and suitable for **educational heat
transfer simulations**.

------------------------------------------------------------------------

# Possible Future Improvements

Potential extensions include:

-   Successive Over‑Relaxation (SOR) for faster convergence
-   Transient heat conduction solver
-   Non‑uniform grids
-   Neumann boundary conditions
-   Interactive GUI visualization

------------------------------------------------------------------------

# Author

Developed as part of **Heat and Mass Transfer numerical simulations**.
