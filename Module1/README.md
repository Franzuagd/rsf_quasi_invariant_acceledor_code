Linear Lattice Module
`linear_lattice.py` builds the complete accelerator ring, calculates its linear optics, performs first-order chromatic correction, and efficiently updates the prepared ring when optimizer parameters change.
The module is designed around two public operations:
```python
prepare_lattice()
update_linear()
```
The normal workflow is to prepare the full ring once and then use `update_linear()` during nonlinear calculations or optimization.
Requirements
Install the required packages with:
```bash
pip install numpy matplotlib
```
NumPy is required by the core module. Matplotlib is only imported when the optional plotting function is used.
Project Files
```text
linear_lattice.py
    Full-ring definition, linear-optics calculations, chromatic correction,
    and efficient update logic.

test_linear_lattice.py
    Optional validation report for the full ring and update behavior.
```
Design
The module follows three rules:
The complete ring is always built and analyzed.
Every lattice-specific definition that a user may need to edit is contained in section `1.1.1`.
After initialization, only calculations affected by edited parameters are repeated.
The magnet representation is a plain Python list:
```python
[name, type, length, angle, K, S, O, M, M5]
```
Repeated occurrences of one magnet family share the same object. Updating a family therefore updates every occurrence of that family in the ring.
User-Editable Section
All manual lattice configuration is located under:
```python
# 1.1.1 USER-EDITABLE FULL-RING DEFINITION
```
This section contains:
the initial parameter dictionary,
the optimizer variable list,
all magnet definitions,
the lattice ordering and number of repetitions,
the full-ring definition,
the linear and chromatic variable categories,
the parameter-to-magnet mapping,
the chromatic correction families and targets,
the default sampling step.
When adapting the module to another lattice, edit section `1.1.1` only. The calculation and update code below that section should not require lattice-specific changes.
Preparing the Full Ring
Import the module and prepare the lattice with:
```python
import linear_lattice as lin

magnets, lattice, data, correction, parameters = lin.prepare_lattice()
```
`prepare_lattice()` performs the complete initialization:
```text
read the parameter dictionary
        ↓
define the unique magnets
        ↓
build the complete ordered ring
        ↓
calculate periodic Twiss and dispersion
        ↓
calculate tunes, chromaticities, emittance, and circumference
        ↓
solve the first-order chromatic correction
```
The current ring contains:
```text
20 repeated cells
2,220 ordered elements
```
There is no cell/ring selection. `prepare_lattice()` always returns the full ring and all ring-level calculations use a repetition factor of one.
Returned Values
Value	Meaning
`magnets`	Unique magnet-family objects
`lattice`	Complete ordered ring
`data`	Linear-optics results and sampled functions
`correction`	Chromatic correction result, or `None` if disabled
`parameters`	Parameter dictionary synchronized with the prepared ring
The chromatic correction has the form:
```python
[
    family1,
    strength1,
    family2,
    strength2,
    corrected_chrom_x,
    corrected_chrom_y,
]
```
Important linear results can be accessed with:
```python
data[lin.CS0]
data[lin.DISP0]
data[lin.LATTICE_M4]
data[lin.LATTICE_M5]
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
`CELL_M4` and `CELL_M5` remain available as compatibility aliases for `LATTICE_M4` and `LATTICE_M5`.
Updating the Prepared Ring
After initialization, use `update_linear()` instead of rebuilding the ring during every optimizer evaluation.
The correct order is:
Change values in the returned `parameters` dictionary.
Record the exact names of the changed parameters.
Call `update_linear()`.
Use the returned `lattice`, `data`, and `parameters`.
For example:
```python
parameters["X1"] = new_X1
parameters["ko1"] = new_ko1

edited_variables = ["X1", "ko1"]

lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    edited_variables,
)

if new_correction is not None:
    correction = new_correction
```
The dictionary contains the new values. `edited_variables` tells the function which values must be copied into the ring and which calculations are no longer valid.
Changing the dictionary without listing the parameter does not update the magnet:
```python
parameters["X1"] = new_X1

# Incorrect: X1 is missing from edited_variables.
lin.update_linear(lattice, data, parameters, [])
```
Listing a parameter before changing its dictionary value simply writes the old value back into the ring.
Recalculation Rules
Changed parameters	Magnet values	Linear optics	Chromatic correction
`energy`, `LSD`, `F`, `X1`–`X9`	Update	Recompute	Recompute
Sextupole parameters	Update	Reuse	Recompute
`ko1`, `ko2`, `ko3`	Update	Reuse	Reuse
The dependency logic is:
```text
linear change
    → update the affected magnets and matrices
    → recompute linear optics
    → recompute chromatic correction

sextupole change
    → update sextupole strengths
    → reuse Twiss and dispersion
    → recompute chromatic correction

thin-multipole change
    → update nonlinear strengths
    → reuse linear optics
    → keep the existing chromatic correction
```
When variables from multiple categories are passed together, the function performs every calculation required by the most influential category.
Understanding `new_correction`
For a nonlinear-only update:
```python
new_correction is None
```
This means that chromatic correction was not recalculated. It does not mean that the corrected sextupole strengths were removed from the ring.
Preserve the previous correction information with:
```python
lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    edited_variables,
)

if new_correction is not None:
    correction = new_correction
```
Updating an Optimizer Candidate
An optimizer can provide a dictionary containing one candidate:
```python
candidate = {
    "X1": new_X1,
    "X2": new_X2,
    "ko1": new_ko1,
}

for name, value in candidate.items():
    parameters[name] = value

lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    list(candidate),
)

if new_correction is not None:
    correction = new_correction
```
The nonlinear calculation should always receive the returned `lattice` and `data`, because linear-variable changes invalidate the previous optics.
Chromatic Correction Families
The default correction families and targets are defined in section `1.1.1`:
```python
CHROMATIC_FAMILY1 = "SF1"
CHROMATIC_FAMILY2 = "SD1"
TARGET_CHROM_X = 0.0
TARGET_CHROM_Y = 0.0
```
When chromatic correction is enabled, the solved correction-family strengths are written into the ring and synchronized with the parameter dictionary through `CORRECTION_PARAMETER_MAP`.
If a correction-family strength is supplied as an edited variable, it may immediately be replaced by the strength required to reach the target chromaticities. To preserve an explicitly supplied value, disable the correction for that update:
```python
lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    ["ksf1", "ksd1"],
    correct_chromatic=False,
)
```
Complete Rebuild
Use `prepare_lattice()` again instead of `update_linear()` when:
starting an independent calculation,
changing the ring ordering,
changing the number of repeated cells,
adding or removing magnets,
changing magnet-family names,
the parameter dictionary and ring may no longer be synchronized.
For example:
```python
p = lin.PARAMETERS.copy()
p["X1"] = new_X1

magnets, lattice, data, correction, parameters = lin.prepare_lattice(
    parameters=p
)
```
`update_linear()` changes supported parameter values. It does not change the ring topology.
Interface With `nonlinear.py`
The intended initialization is:
```python
import linear_lattice as lin

magnets, lattice, data, correction, parameters = lin.prepare_lattice()
```
The nonlinear module should use the returned `lattice` and `data` rather than reconstructing the accelerator internally.
During optimization:
```text
optimizer candidate
        ↓
update parameters
        ↓
update_linear()
        ↓
updated full ring + valid linear data
        ↓
nonlinear calculation
```
This keeps the linear and nonlinear models synchronized while avoiding unnecessary optics calculations.
Optional Validation and Plotting
The final section of `linear_lattice.py` contains optional utilities for:
consistency checks,
terminal summaries,
plotting Twiss functions and dispersion,
standalone execution.
Run the module directly with:
```bash
python linear_lattice.py
```
Run the detailed report with:
```bash
python test_linear_lattice.py
```
The validation checks:
symplecticity,
periodic Twiss closure,
periodic dispersion closure,
Courant–Snyder identities,
total ring bending,
nonlinear-only updates,
sextupole-only updates,
linear-optics updates.
The current full-ring results agree with the independent OPA reference values.
Recommended Workflow
```python
import linear_lattice as lin

# Prepare the full ring once.
magnets, lattice, data, correction, parameters = lin.prepare_lattice()

# Update one optimizer candidate.
candidate = {
    "X1": new_X1,
    "ko1": new_ko1,
}

parameters.update(candidate)

lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    candidate.keys(),
)

if new_correction is not None:
    correction = new_correction

# Use the returned lattice and data in the nonlinear calculation.
```
In short:
```text
prepare once
    → edit the parameter dictionary
    → call update_linear()
    → use the returned ring and linear data
```
