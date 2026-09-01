# NOTE
This is README file is generated using chagpt for the moment. Code is still not proving to be working.

# Polynomial Quasi-Invariant Accelerator Optimization

This project implements a modular workflow for studying and optimizing nonlinear accelerator lattices using polynomial quasi-invariants.

The code is organized into three main stages:

```text
linear_lattice.py
        ↓
    nonlinear.py
        ↓
   optimization.py
        ↓
       run.py
```

The project separates linear optics, nonlinear dynamics, and optimization so that each stage can be tested independently.

---

# Project Structure

```text
project/
│
├── linear_lattice.py
├── nonlinear.py
├── optimization.py
├── run.py
│
├── test_linear_lattice.py
└── README.md
```

## `linear_lattice.py`

Builds the accelerator lattice and computes the linear optics.

It provides:

```text
lattice definition
        ↓
transfer matrices
        ↓
periodic Twiss functions
        ↓
dispersion
        ↓
tunes
        ↓
chromaticities
        ↓
natural emittance
        ↓
chromatic correction
```

The main entry point is:

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

The important outputs are:

```text
cell
    ordered accelerator cell

data
    linear optics information

correction
    corrected sextupole strengths and chromaticities

parameters
    parameters used to build the lattice
```

The magnet representation used throughout the project is:

```python
[name, type, length, angle, K, S, O, M, M5]
```

---

# `nonlinear.py`

Uses the lattice and linear optics produced by `linear_lattice.py`.

The general workflow is:

```text
linear lattice
      │
      ▼
polynomial basis
      │
      ▼
Hamiltonian representation
      │
      ▼
nonlinear element maps
      │
      ▼
one-cell nonlinear transfer matrix
      │
      ▼
Courant-Snyder quadratic invariant
      │
      ▼
higher-order quasi-invariant
```

Typical use:

```python
import nonlinear as nl

Ix, Iy, result, transfer = nl.invariant_vectors(
    cell,
    data,
    order=8,
    delta_order=1,
)
```

The main outputs are:

```text
Ix
    horizontal quasi-invariant vector

Iy
    vertical quasi-invariant vector

transfer
    nonlinear one-cell transfer matrix

result
    numerical information from the invariant construction
```

For visualization or inspection, the invariant vectors can also be converted into symbolic polynomials or plotted.

---

# `optimization.py`

Connects the linear and nonlinear modules to an optimizer.

For every candidate parameter vector, the workflow is:

```text
candidate parameters
        │
        ▼
build new lattice
        │
        ▼
compute linear optics
        │
        ▼
correct chromaticity
        │
        ▼
build nonlinear transfer map
        │
        ▼
construct quasi-invariant
        │
        ▼
evaluate objective
```

The optimizer therefore always evaluates a complete, self-consistent lattice candidate.

The optimization uses a hybrid strategy:

```text
CMA-ES
   ↓
Powell
```

CMA-ES performs the exploratory search and Powell refines the best candidate.

---

# `run.py`

`run.py` is the main user-facing script for long optimization runs.

Run it with:

```bash
python run.py
```

A typical production configuration uses:

```text
polynomial order: m = 8
momentum order:   d = 1
```

The run script controls:

```text
optimization variables
CMA-ES settings
Powell settings
linear constraints
runtime
output folders
plots
```

---

# Complete Workflow

The complete project flow is:

```text
                PARAMETERS
                    │
                    ▼
        ┌─────────────────────┐
        │  linear_lattice.py  │
        └─────────────────────┘
                    │
                    │ cell + linear optics
                    ▼
          ┌──────────────────┐
          │   nonlinear.py   │
          └──────────────────┘
                    │
                    │ quasi-invariant
                    ▼
        ┌─────────────────────┐
        │   optimization.py   │
        └─────────────────────┘
                    │
                    │ objective value
                    ▼
               optimizer
                    │
                    ▼
             optimized lattice
```

For each optimization candidate:

```text
parameter vector
      ↓
parameter dictionary
      ↓
linear lattice
      ↓
periodic optics
      ↓
chromatic correction
      ↓
nonlinear map
      ↓
quasi-invariant
      ↓
objective value
```

---

# Basic Usage

## Linear lattice only

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()

print(data[lin.TUNE_X])
print(data[lin.TUNE_Y])
print(data[lin.EMITTANCE])
```

---

## Linear + nonlinear calculation

```python
import linear_lattice as lin
import nonlinear as nl

magnets, cell, data, correction, parameters = lin.prepare_lattice()

Ix, Iy, result, transfer = nl.invariant_vectors(
    cell,
    data,
    order=8,
    delta_order=1,
)
```

---

## Optimization

```python
import numpy as np
import optimization as opt

vary = list(opt.VARY)

v0 = np.array(
    [
        opt.BASE_PARAMETERS[name]
        for name in vary
    ],
    dtype=float,
)

vfinal, ffinal, es, res, history = opt.hybrid_optimize(
    v0
)
```

---

# Testing

The linear lattice can be tested independently using:

```bash
python test_linear_lattice.py
```

The test checks quantities such as:

```text
symplecticity
Twiss closure
dispersion closure
Courant-Snyder identities
ring bending
tunes
chromaticities
emittance
```

The nonlinear module can also be tested independently:

```bash
python nonlinear.py
```

Using smaller polynomial orders is recommended for quick development tests.

---

# Validation

The project is designed so each stage can be validated separately.

The linear lattice can be compared against an external accelerator optics program using:

```text
beta functions
alpha functions
dispersion
tunes
circumference
chromaticities
emittance
```

The nonlinear and optimization results should then be validated independently with particle tracking and dynamic-aperture calculations.

---

# Requirements

Install the required packages with:

```bash
pip install numpy scipy sympy matplotlib cma
```

Main dependencies:

```text
NumPy
SciPy
SymPy
Matplotlib
CMA-ES
```

---

# Summary

The project can be viewed as:

```text
accelerator lattice
        ↓
linear optics
        ↓
nonlinear polynomial map
        ↓
quasi-invariant
        ↓
optimization objective
        ↓
optimized accelerator parameters
```

The main design goal is to keep each part of the calculation modular, testable, and reusable.
