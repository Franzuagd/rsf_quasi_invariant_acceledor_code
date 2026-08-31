#IMPORTANT
THIS README WAS CREATED USING CHAGPT.
# Linear Lattice Module 

`linear_lattice.py` provides the linear accelerator-lattice model used by the project.

The module includes:

* magnet and lattice data structures,
* linear transfer matrices,
* periodic Courant–Snyder/Twiss functions,
* horizontal dispersion,
* tunes and natural chromaticities,
* natural horizontal emittance,
* first-order chromaticity correction,
* lattice consistency checks,
* optional plotting and summary helpers.

The module is designed to be imported by later nonlinear and optimization code.

---

## Requirements

The module requires:

```bash
pip install numpy
```

`matplotlib` is only required when using `plot_linear_functions()`:

```bash
pip install matplotlib
```

---

# Basic Usage

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

This is the recommended entry point.

It returns:

```text
magnets     unique magnet definitions
cell        ordered accelerator cell
data        linear optics data
correction  chromaticity-correction result, or None
parameters  copy of the parameters used to build the lattice
```

For example:

```python
print(data[lin.TUNE_X])
print(data[lin.TUNE_Y])
print(data[lin.EMITTANCE])
print(data[lin.CIRCUMFERENCE])
```

---

# 1. Magnet Data Structure

Every magnet is stored as a Python list:

```python
[name, type, length, angle, K, S, O, M, M5]
```

The corresponding indices are:

```python
lin.NAME
lin.TYPE
lin.LENGTH
lin.ANGLE
lin.K
lin.S
lin.O
lin.M
lin.M5
```

Example:

```python
qf1 = lin.get_magnet(cell, "QF1")

print(qf1[lin.NAME])
print(qf1[lin.LENGTH])
print(qf1[lin.K])
```

## Fields

| Field    | Meaning                                                                 |
| -------- | ----------------------------------------------------------------------- |
| `NAME`   | Magnet or family name                                                   |
| `TYPE`   | `"drift"`, `"quadrupole"`, `"sextupole"`, `"bending"`, or `"multipole"` |
| `LENGTH` | Physical element length                                                 |
| `ANGLE`  | Total bending angle in degrees                                          |
| `K`      | Quadrupole strength                                                     |
| `S`      | Sextupole strength                                                      |
| `O`      | Thin higher-multipole kick strength                                     |
| `M`      | Cached 4×4 linear transfer matrix                                       |
| `M5`     | Cached 5×5 map used for dispersion                                      |

Thin `multipole` elements have zero physical length and no linear effect.

---

# 2. Parameter Dictionary

The default lattice parameters are stored in:

```python
lin.PARAMETERS
```

Example:

```python
p = lin.PARAMETERS.copy()

print(p["energy"])
print(p["X1"])
print(p["ks1"])
```

The dictionary contains the beam energy, lattice geometry parameters, quadrupole strengths, sextupole strengths, and thin multipole strengths.

To build the lattice using modified parameters:

```python
p = lin.PARAMETERS.copy()

p["X1"] = 3.7
p["X2"] = -4.3

magnets = lin.define_magnets(p)
cell = lin.build_lattice(lin.CELL_NAMES, magnets)

data = lin.linear_optics(
    cell,
    energy=p["energy"],
    repetitions=lin.N_CELLS,
)
```

---

# 3. Lattice Structure

The ordered cell is defined by:

```python
lin.CELL_NAMES
```

The number of repeated cells is:

```python
lin.N_CELLS
```

An explicit full-ring name list is also available:

```python
lin.RING_NAMES
```

Build the cell with:

```python
magnets = lin.define_magnets(lin.PARAMETERS)
cell = lin.build_lattice(lin.CELL_NAMES, magnets)
```

`build_lattice()` preserves the exact order of `CELL_NAMES`.

Repeated occurrences of the same magnet family reference the same magnet object.

Therefore, changing one family updates every occurrence of that family in the cell.

---

# 4. Updating Magnet Strengths

## Quadrupoles

Use:

```python
lin.set_quadrupole_strength(
    cell,
    "QF1",
    3.8,
)
```

This:

1. changes the quadrupole strength `K`,
2. recomputes the corresponding 4×4 matrix,
3. recomputes the corresponding 5×5 matrix.

After changing a quadrupole, recompute the optics:

```python
data = lin.linear_optics(
    cell,
    energy=parameters["energy"],
    repetitions=lin.N_CELLS,
)
```

## Sextupoles

Use:

```python
lin.set_sextupole_strength(
    cell,
    "S1",
    -30.0,
)
```

The sextupole linear map is treated as a drift, so changing `S` does not require recomputing `M` or `M5`.

## Manual Changes

If a linear parameter is modified directly:

```python
elem = lin.get_magnet(cell, "QF1")
elem[lin.K] = 3.8
```

refresh the matrices with:

```python
lin.refresh_matrices(cell)
```

When possible, prefer the provided update functions.

---

# 5. Transfer-Matrix Functions

## `plane_matrix(L, k)`

Returns the 2×2 transfer matrix for constant focusing:

```python
M = lin.plane_matrix(L, k)
```

## `drift4(L)` and `drift5(L)`

```python
M4 = lin.drift4(L)
M5 = lin.drift5(L)
```

## `quadrupole4(L, k)` and `quadrupole5(L, k)`

```python
M4 = lin.quadrupole4(L, k)
M5 = lin.quadrupole5(L, k)
```

## `bending4(L, k, h)` and `bending5(L, k, h)`

```python
M4 = lin.bending4(L, k, h)
M5 = lin.bending5(L, k, h)
```

where:

```text
h = 1 / rho
```

is the curvature.

## `transfer4(elem, L=None)`

Returns the appropriate 4×4 map for a lattice element:

```python
M4 = lin.transfer4(elem)
```

A partial element can also be propagated:

```python
M4 = lin.transfer4(elem, 0.01)
```

## `transfer5(elem, L=None)`

Same idea, but using the 5×5 representation:

```python
M5 = lin.transfer5(elem)
```

## `compute_matrices(elem)`

Recomputes and stores the full-element matrices:

```python
lin.compute_matrices(elem)
```

## `lattice_matrix(lattice)`

Returns the complete transfer matrices of an ordered lattice:

```python
M4, M5 = lin.lattice_matrix(cell)
```

---

# 6. Periodic Twiss and Dispersion

Use:

```python
cs0, disp0, M4, M5 = lin.periodic_twiss_and_dispersion(cell)
```

`cs0` contains:

```text
[
    beta_x,
    alpha_x,
    gamma_x,
    beta_y,
    alpha_y,
    gamma_y,
]
```

`disp0` contains:

```text
[
    Dx,
    Dpx,
    0,
    0,
    1,
]
```

where `Dx` is the periodic horizontal dispersion.

---

# 7. Complete Linear Optics

Use:

```python
data = lin.linear_optics(
    cell,
    energy=3.0,
    repetitions=lin.N_CELLS,
    step=0.01,
)
```

The result is accessed using named index constants.

| Constant        | Meaning                                 |
| --------------- | --------------------------------------- |
| `CS0`           | Periodic Courant–Snyder parameters      |
| `DISP0`         | Periodic dispersion                     |
| `CELL_M4`       | One-cell 4×4 matrix                     |
| `CELL_M5`       | One-cell 5×5 matrix                     |
| `TUNE_X`        | Full-ring horizontal tune               |
| `TUNE_Y`        | Full-ring vertical tune                 |
| `CHROM_X`       | Natural horizontal chromaticity         |
| `CHROM_Y`       | Natural vertical chromaticity           |
| `EMITTANCE`     | Natural horizontal emittance            |
| `CIRCUMFERENCE` | Full-ring circumference                 |
| `S_VALUES`      | Longitudinal positions through one cell |
| `CS_VALUES`     | Twiss functions along one cell          |
| `DISP_VALUES`   | Dispersion along one cell               |

Example:

```python
print("Qx =", data[lin.TUNE_X])
print("Qy =", data[lin.TUNE_Y])

print("Chromaticities:")
print(data[lin.CHROM_X])
print(data[lin.CHROM_Y])

print("Emittance =", data[lin.EMITTANCE])
print("Circumference =", data[lin.CIRCUMFERENCE])
```

---

# 8. Linear Functions Along the Cell

The lower-level propagation function can also be used directly:

```python
cs0, disp0, _, _ = lin.periodic_twiss_and_dispersion(cell)

result = lin.propagate_linear_functions(
    cell,
    cs0,
    disp0,
    step=0.01,
)
```

It returns:

```text
[
    mux,
    muy,
    s_values,
    cs_values,
    disp_values,
    radiation,
    chromx_integral,
    chromy_integral,
]
```

For normal use, `linear_optics()` is preferred because it converts these quantities into full-ring tunes, chromaticities, emittance, and circumference.

---

# 9. Chromaticity Correction

The module supports first-order chromaticity correction using two sextupole families.

The default families are:

```text
SF1
SD1
```

The correction uses the effective focusing:

```text
K_eff = K - 2 Dx S
```

## Recommended Usage

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice(
    correct_chromatic=True,
    target_chrom_x=0.0,
    target_chrom_y=0.0,
)
```

The returned `correction` contains:

```text
[
    family1,
    strength1,
    family2,
    strength2,
    corrected_chrom_x,
    corrected_chrom_y,
]
```

Example:

```python
print(correction)
```

The corrected sextupole strengths are also stored directly in the shared magnet objects.

## Disable Chromatic Correction

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice(
    correct_chromatic=False
)
```

Then:

```python
correction is None
```

## Different Targets

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice(
    target_chrom_x=1.0,
    target_chrom_y=1.0,
)
```

## Available Correction Functions

```python
lin.chromaticity_with_sextupoles(...)
lin.correct_chromaticity(...)
lin.correct_chromaticity_from_data(...)
```

`correct_chromaticity_from_data()` is useful when the linear optics has already been computed because it reuses the stored Twiss and dispersion arrays.

---

# 10. Finding Magnets

Use:

```python
qf1 = lin.get_magnet(cell, "QF1")
```

Example:

```python
print(qf1[lin.K])
```

If the magnet does not exist, `get_magnet()` raises a `KeyError`.

To retrieve each shared magnet object only once:

```python
unique = lin.unique_magnets(cell)
```

---

# 11. Consistency Checks

Use:

```python
checks = lin.check_linear_lattice(
    cell,
    data,
)
```

Available check indices are:

```python
lin.CHECK_SYMPLECTIC
lin.CHECK_TWISS_CLOSURE
lin.CHECK_DISPERSION_CLOSURE
lin.CHECK_CS_X
lin.CHECK_CS_Y
lin.CHECK_RING_BENDING_DEG
```

Example:

```python
print(
    checks[lin.CHECK_SYMPLECTIC]
)

print(
    checks[lin.CHECK_TWISS_CLOSURE]
)

print(
    checks[lin.CHECK_RING_BENDING_DEG]
)
```

The checks include:

* 4D symplecticity,
* periodic Twiss closure,
* periodic dispersion closure,
* Courant–Snyder identities,
* total ring bending angle.

---

# 12. Display Helpers

Print the main quantities:

```python
lin.print_linear_summary(
    data,
    parameters["energy"],
    correction,
)
```

Plot the linear functions:

```python
lin.plot_linear_functions(data)
```

The plot includes:

```text
beta_x(s)
beta_y(s)
Dx(s)
```

through one cell.

---

# 13. Recommended Update Workflows

## Change One Quadrupole

```python
lin.set_quadrupole_strength(
    cell,
    "QF1",
    new_k,
)

data = lin.linear_optics(
    cell,
    energy=parameters["energy"],
    repetitions=lin.N_CELLS,
)
```

## Change One Sextupole

```python
lin.set_sextupole_strength(
    cell,
    "S1",
    new_s,
)
```

If corrected chromaticities are required afterward, run the chromaticity correction again.

## Change Several Parameters

```python
p = lin.PARAMETERS.copy()

p["X1"] = new_x1
p["X2"] = new_x2
p["ks1"] = new_ks1

magnets = lin.define_magnets(p)

cell = lin.build_lattice(
    lin.CELL_NAMES,
    magnets,
)

data = lin.linear_optics(
    cell,
    energy=p["energy"],
    repetitions=lin.N_CELLS,
)
```

## Rebuild Everything

The simplest option when several parameters change is:

```python
p = lin.PARAMETERS.copy()

p["X1"] = new_x1
p["X2"] = new_x2

magnets, cell, data, correction, parameters = lin.prepare_lattice(
    parameters=p
)
```

---

# 14. Interface With the Nonlinear Module

The recommended interface is:

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

The nonlinear module can then use:

```python
cell
data
```

Magnet information is available using:

```python
elem[lin.TYPE]
elem[lin.LENGTH]
elem[lin.ANGLE]
elem[lin.K]
elem[lin.S]
elem[lin.O]
```

Linear optics information is available using:

```python
data[lin.CS0]
data[lin.DISP0]
data[lin.CELL_M4]
data[lin.CELL_M5]
data[lin.CS_VALUES]
data[lin.DISP_VALUES]
```

---

# 15. Standalone Use

The module can be executed directly:

```bash
python linear_lattice.py
```

This:

* builds the default lattice,
* computes the linear optics,
* performs chromaticity correction,
* prints the main accelerator parameters,
* runs the consistency checks.

Importing the module does not execute the standalone test:

```python
import linear_lattice
```

because the executable section is protected by:

```python
if __name__ == "__main__":
    main()
```

---

# Main Public Functions

For quick reference:

```python
# Lattice construction
magnet(...)
define_magnets(...)
build_lattice(...)
get_magnet(...)

# Magnet updates
set_quadrupole_strength(...)
set_sextupole_strength(...)
refresh_matrices(...)

# Transfer matrices
plane_matrix(...)
drift4(...)
drift5(...)
quadrupole4(...)
quadrupole5(...)
bending4(...)
bending5(...)
transfer4(...)
transfer5(...)
lattice_matrix(...)

# Linear optics
periodic_twiss_and_dispersion(...)
propagate_linear_functions(...)
linear_optics(...)

# Chromaticity
chromaticity_with_sextupoles(...)
correct_chromaticity(...)
correct_chromaticity_from_data(...)

# Main frontend
prepare_lattice(...)

# Verification
check_linear_lattice(...)

# Display
print_linear_summary(...)
plot_linear_functions(...)
```

For most users, the main workflow is simply:

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```
