# Linear Lattice Module

`linear_lattice.py` contains the linear accelerator model used by the project.

It provides:

* lattice construction,
* linear transfer matrices,
* periodic Twiss functions,
* horizontal dispersion,
* tunes,
* natural chromaticities,
* natural horizontal emittance,
* first-order chromaticity correction,
* linear consistency checks.

The module uses a plain Python-list representation for magnets:

```python
[name, type, length, angle, K, S, O, M, M5]
```

---

# First Run

Install the required packages:

```bash
pip install numpy matplotlib
```

Run the module directly with:

```bash
python linear_lattice.py
```

or run the dedicated test:

```bash
python test_linear_lattice.py
```

The test builds the lattice, calculates the linear optics, performs the chromatic correction, checks the lattice consistency, and saves a text report.

For use from another module:

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

This is the recommended entry point.

---

# Main Output

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

returns:

```text
magnets
    unique magnet definitions

cell
    ordered accelerator cell

data
    calculated linear optics

correction
    chromaticity-correction information

parameters
    parameter dictionary used to construct the lattice
```

Typical values can be accessed with:

```python
print(data[lin.TUNE_X])
print(data[lin.TUNE_Y])
print(data[lin.CHROM_X])
print(data[lin.CHROM_Y])
print(data[lin.EMITTANCE])
print(data[lin.CIRCUMFERENCE])
```

---

# User-Level Configuration

The main lattice parameters are stored in:

```python
lin.PARAMETERS
```

For example:

```python
PARAMETERS = {
    "energy": 3.0,
    "LSD": 0.1,
    "F": 0.8,

    "X1": ...,
    "X2": ...,
    "X3": ...,

    "ks1": ...,
    "ks2": ...,

    "ko1": ...,
    "ko2": ...,
    "ko3": ...,
}
```

These values describe the accelerator and may be changed by the user.

The names selected for later optimization are stored in:

```python
lin.VARY
```

At the moment, `VARY` should be treated only as the current list of parameters intended to be varied. The actual optimization workflow is handled outside this module.

Other important global configuration values are:

```python
lin.CELL_NAMES
lin.N_CELLS
lin.RING_NAMES
```

These define the lattice ordering and number of repeated cells.

---

# Internal Constants

The magnet-list positions are defined by:

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

These are indices used to access the magnet data structure.

For example:

```python
qf1 = lin.get_magnet(cell, "QF1")

print(qf1[lin.K])
print(qf1[lin.LENGTH])
```

`M` and `M5` are internal cached linear maps:

```text
M   4×4 transverse linear map
M5  5×5 map used for dispersion
```

They should normally be updated through the module functions rather than edited directly.

---

# Magnet Data Structure

Each magnet has the form:

```python
[
    name,
    type,
    length,
    angle,
    K,
    S,
    O,
    M,
    M5,
]
```

with:

| Field    | Meaning                       |
| -------- | ----------------------------- |
| `NAME`   | Magnet/family name            |
| `TYPE`   | Magnet type                   |
| `LENGTH` | Physical length               |
| `ANGLE`  | Total bending angle           |
| `K`      | Quadrupole strength           |
| `S`      | Sextupole strength            |
| `O`      | Thin nonlinear multipole kick |
| `M`      | Cached 4×4 linear map         |
| `M5`     | Cached 5×5 dispersion map     |

Supported types include:

```text
drift
quadrupole
sextupole
bending
multipole
```

Repeated occurrences of the same family in the cell share the same magnet object.

Therefore, changing one family changes every occurrence of that family.

---

# Building the Lattice

The simplest way to build everything is:

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

The workflow is:

```text
PARAMETERS
    ↓
define magnets
    ↓
build ordered cell
    ↓
compute periodic linear optics
    ↓
compute Twiss and dispersion along the cell
    ↓
compute tunes, chromaticities and emittance
    ↓
solve chromatic correction
```

The cell and optics are then ready to be passed to the nonlinear module.

---

# Using Different Parameters

The safest general method is to create a new parameter dictionary:

```python
p = lin.PARAMETERS.copy()

p["X1"] = new_value
p["X2"] = new_value

magnets, cell, data, correction, parameters = lin.prepare_lattice(
    parameters=p
)
```

`prepare_lattice()` builds a fresh lattice from the supplied parameters.

Therefore, when this function is called again, the module currently recomputes:

```text
magnet definitions
cell
linear matrices
periodic Twiss
dispersion
radiation integrals
emittance
natural chromaticity
chromatic correction
```

This is the reliable general update path.

---

# Updating an Existing Cell

The module also provides direct update functions when a complete rebuild is not desired.

## Quadrupole

```python
lin.set_quadrupole_strength(
    cell,
    "QF1",
    new_K,
)
```

Changing a quadrupole changes its linear transfer matrix.

The function therefore updates:

```text
K
M
M5
```

After changing a quadrupole, the stored `data` is no longer valid.

Recompute the linear optics:

```python
data = lin.linear_optics(
    cell,
    energy=parameters["energy"],
    repetitions=lin.N_CELLS,
)
```

---

# Updating Sextupoles

A sextupole can be changed with:

```python
lin.set_sextupole_strength(
    cell,
    "S1",
    new_S,
)
```

In the current linear model, the on-momentum sextupole linear map is treated as a drift.

Therefore changing `S` does **not** change:

```text
M
M5
Twiss functions
dispersion
tunes
emittance
```

So the complete linear optics does not need to be recalculated only because a sextupole strength changed.

However, the first-order chromaticity **does** change.

---

# Chromaticity Correction

The default chromatic correction uses:

```text
SF1
SD1
```

as the two correction families.

The normal full calculation is:

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice(
    correct_chromatic=True,
    target_chrom_x=0.0,
    target_chrom_y=0.0,
)
```

This solves the required `SF1` and `SD1` strengths automatically.

The corrected strengths are written directly into the shared magnet objects in `cell`.

---

# Updating Chromaticity After Sextupole Changes

If the linear optics has **not changed**, it is unnecessary to propagate the Twiss functions and dispersion again.

For example:

```python
lin.set_sextupole_strength(
    cell,
    "S1",
    new_S,
)
```

can be followed by:

```python
S1, S2, chrom_x, chrom_y = lin.correct_chromaticity_from_data(
    cell,
    data,
    family1="SF1",
    family2="SD1",
    target_x=0.0,
    target_y=0.0,
    repetitions=lin.N_CELLS,
)
```

This reuses:

```python
data[lin.CS_VALUES]
data[lin.DISP_VALUES]
```

that were already calculated.

The workflow is therefore:

```text
existing cell + existing linear data
              ↓
change sextupole strengths
              ↓
reuse Twiss + dispersion
              ↓
solve new SF1 / SD1 strengths
```

This avoids an unnecessary second linear-optics calculation.

---

# When Must Chromaticity Correction Be Recomputed?

If a quantity that changes the linear optics is modified, for example:

```text
quadrupole strength
bend
magnet length
lattice geometry
```

then:

```text
Twiss changes
dispersion changes
        ↓
sextupole chromatic response changes
        ↓
SF1 / SD1 must be solved again
```

In this situation the easiest and safest method is simply:

```python
magnets, cell, data, correction, parameters = lin.prepare_lattice(
    parameters=p
)
```

---

# If Only Sextupoles Change

The efficient workflow is:

```text
keep existing linear data
        ↓
change sextupole S
        ↓
recompute chromatic correction
        ↓
pass corrected cell to nonlinear.py
```

The linear Twiss and dispersion do not need to be recalculated.

---

# If Only Thin Multipoles Change

The `O` value represents the thin nonlinear multipole strength.

Thin multipoles have no linear effect in this module.

Therefore changing only `O` does not require recomputing:

```text
Twiss
dispersion
tunes
emittance
first-order chromaticity
```

The updated `cell` can be passed directly to the nonlinear calculation.

The current module does not provide a dedicated `set_multipole_strength()` helper, so a thin multipole can currently be updated with:

```python
elem = lin.get_magnet(cell, "O1")
elem[lin.O] = new_O
```

No `M` or `M5` refresh is required because a thin multipole has an identity linear map.

---

# Recalculation Rules

For the current module:

| Changed quantity   | Linear optics | Chromatic correction |
| ------------------ | ------------: | -------------------: |
| Quadrupole `K`     |     Recompute |            Recompute |
| Bending / geometry |     Recompute |            Recompute |
| Sextupole `S`      |         Reuse |            Recompute |
| Thin multipole `O` |         Reuse |                Reuse |

In short:

```text
linear change
    → recompute linear optics
    → recompute chromatic correction

sextupole change
    → reuse linear optics
    → recompute chromatic correction

thin multipole change
    → reuse linear optics
    → reuse chromatic correction
```

---

# Current `prepare_lattice()` Behavior

It is important to distinguish what is **possible** from what `prepare_lattice()` currently does.

Calling:

```python
lin.prepare_lattice(parameters=p)
```

always builds a fresh lattice and runs the linear-optics calculation.

It does not currently inspect which parameter changed and selectively skip calculations.

Therefore:

```text
prepare_lattice()
    = safe complete rebuild
```

For one-off calculations this is normally the easiest method.

For performance-sensitive optimization, the more selective update methods described above can be used when the linear lattice is known to remain unchanged.

---

# Linear Data

The main calculated data is available through:

```python
data[lin.CS0]
data[lin.DISP0]

data[lin.CELL_M4]
data[lin.CELL_M5]

data[lin.TUNE_X]
data[lin.TUNE_Y]

data[lin.CHROM_X]
data[lin.CHROM_Y]

data[lin.EMITTANCE]
data[lin.CIRCUMFERENCE]

data[lin.S_VALUES]
data[lin.CS_VALUES]
data[lin.DISP_VALUES]
```

`CS_VALUES` contains:

```text
beta_x
alpha_x
gamma_x
beta_y
alpha_y
gamma_y
```

along the cell.

`DISP_VALUES` contains the propagated dispersion data.

---

# Interface With `nonlinear.py`

The intended interface is:

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

followed by passing:

```python
cell
data
```

to the nonlinear module.

`linear_lattice.py` is responsible for preparing the accelerator and its linear optics.

`nonlinear.py` should use the prepared lattice rather than reconstructing the linear accelerator internally.

If only nonlinear magnet strengths change, the existing linear `data` can be reused whenever those changes do not affect the linear optics.

---

# Validation

Run:

```bash
python test_linear_lattice.py
```

The validation report checks:

```text
symplecticity
periodic Twiss closure
periodic dispersion closure
Courant-Snyder identities
total ring bending
tunes
chromaticity
emittance
```

The main linear quantities can also be compared against an independent accelerator-optics code.

---

# Typical Use

For most one-off calculations:

```python
import linear_lattice as lin

magnets, cell, data, correction, parameters = lin.prepare_lattice()
```

For a modified linear lattice:

```python
p = lin.PARAMETERS.copy()

p["X1"] = new_X1

magnets, cell, data, correction, parameters = lin.prepare_lattice(
    parameters=p
)
```

For a sextupole-only change:

```python
lin.set_sextupole_strength(
    cell,
    "S1",
    new_S,
)

S1, S2, chrom_x, chrom_y = lin.correct_chromaticity_from_data(
    cell,
    data,
)
```

For a thin-multipole-only change:

```python
O1 = lin.get_magnet(cell, "O1")
O1[lin.O] = new_O1
```

and the existing linear optics can be reused.

---

# Summary

The module has two possible usage styles.

## Safe full rebuild

```python
lin.prepare_lattice(parameters=p)
```

Use this whenever linear parameters may have changed.

## Selective update

```text
quadrupole change
    → update K
    → recompute linear optics
    → redo chromatic correction

sextupole change
    → update S
    → reuse linear optics
    → redo chromatic correction

thin multipole change
    → update O
    → reuse linear optics
```

The current module already supports these selective operations, but `prepare_lattice()` itself always performs the complete rebuild.

