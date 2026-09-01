linear_lattice(4).py
Python
okay implement my two request by touching as less as posible, do not add or remove coments

This chat now uses Work.Undo

linear_lattice(5).py
Python
Sorry, I pasted the wrong firle, do it for this one
test_linear_lattice(1).py
Python
edit this accordingly
Pasted markdown(1).md
File
now you can addapt this, also explain how to use update properly: 
linear_lattice(6).py
Python
just a little change to update: my goal is that changing the lattice got to be simple, so linear_varaibles, chromatic_varaibles, and parameter_map should be at the biggining once we defined the lattice for easy access, maybe you can make them global parameter or as them as inut in update_linear. In generar, update_linear should be extremely generic so that it works for any lattice, and all lattice dependt varaibles hsould be given as input and specify by the user.

Refactored update_linear() to be lattice-independent.

The following configuration is now located immediately after the lattice definition:

LINEAR_VARIABLES
CHROMATIC_VARIABLES
PARAMETER_MAP
CORRECTION_PARAMETER_MAP
NONLINEAR_VARIABLES
ENERGY_PARAMETER

NONLINEAR_VARIABLES is derived automatically from variables in PARAMETER_MAP that are neither linear nor chromatic.

The current lattice still uses the simple call:

lattice, data, correction, parameters = update_linear(
    lattice,
    data,
    parameters,
    edited_variables,
)

For another lattice, provide its configuration:

lattice, data, correction, parameters = update_linear(
    lattice,
    data,
    parameters,
    edited_variables,
    repetitions=my_repetitions,
    linear_variables=my_linear_variables,
    chromatic_variables=my_chromatic_variables,
    parameter_map=my_parameter_map,
    magnet_builder=my_define_magnets,
    correction_parameter_map=my_correction_map,
    energy_parameter="energy",
)

The mapping format is:

PARAMETER_MAP = {
    "X1": [("QF1", K)],
    "LSD": [
        ("D1", LENGTH),
        ("SE1", LENGTH),
    ],
    "ko1": [("O1", O)],
}

An arbitrary custom-lattice configuration and both cell/ring reports passed testing.

Download 
linear_lattice(6).py

linear_lattice(6).py
Code
two more things: 1. remove the cell ring choice, we are gonna work with the full ring only for simplicity, in the comments make suer everything the user is going to edit is marked in 1.1.1. and after taht there should not be anything the user should edit: I mean lattice definition, magnets definitions, parameters, etc. Reorganized everything neccesary at first, then at the end craft a new section call optional, where you put plotting, testing, and anything that won't likely be used by non_linear or optimized. REmember the main fundamental goal is to load linear data, and be able to updated efficiently.

Reorganized the module around full-ring use only.

Key changes:

prepare_lattice() no longer accepts scope; it always builds the 2,220-element ring.
All editable definitions are contained in 1.1.1, including:
Parameters and optimizer variables
Magnet definitions
Lattice ordering
Update dependency sets and mappings
Chromatic families and targets
Sampling step
Core order is now: matrices → construction → optics → chromaticity → preparation → efficient update.
Testing, reporting, plotting, and main() are grouped under OPTIONAL.
Preserved CELL_M4 and CELL_M5 as compatibility aliases.
Updated the external test for full-ring-only operation.

Normal usage is now:

magnets, lattice, data, correction, parameters = lin.prepare_lattice()

parameters["X1"] = new_X1

lattice, data, correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    ["X1"],
)

Full-ring optics, consistency checks, and all three update categories passed.

linear_lattice(6).py
test_linear_lattice(1).py
linear_lattice(6).py
Code
test_linear_lattice(1).py
Code
now craft the readme file for this new lienar_lattice

Created the README for the full-ring linear_lattice.py.

It covers:

Full-ring-only architecture
User-editable 1.1.1 section
prepare_lattice() workflow and outputs
Correct update_linear() usage
Recalculation dependency rules
Optimizer and nonlinear-module integration
Chromatic correction behavior
Optional testing and plotting
Complete rebuild conditions

All 16 Python examples passed syntax validation.

Download 
README.md

README.md
Document

Library
/
README.md
Linear Lattice Module

linear_lattice.py builds the complete accelerator ring, calculates its linear optics, performs first-order chromatic correction, and efficiently updates the prepared ring when optimizer parameters change.

The module is designed around two public operations:

prepare_lattice()
update_linear()

The normal workflow is to prepare the full ring once and then use update_linear() during nonlinear calculations or optimization.

Requirements

Install the required packages with:

pip install numpy matplotlib

NumPy is required by the core module. Matplotlib is only imported when the optional plotting function is used.

Project Files
linear_lattice.py
    Full-ring definition, linear-optics calculations, chromatic correction,
    and efficient update logic.

test_linear_lattice.py
    Optional validation report for the full ring and update behavior.
Design

The module follows three rules:

The complete ring is always built and analyzed.
Every lattice-specific definition that a user may need to edit is contained in section 1.1.1.
After initialization, only calculations affected by edited parameters are repeated.

The magnet representation is a plain Python list:

[name, type, length, angle, K, S, O, M, M5]

Repeated occurrences of one magnet family share the same object. Updating a family therefore updates every occurrence of that family in the ring.

User-Editable Section

All manual lattice configuration is located under:

# 1.1.1 USER-EDITABLE FULL-RING DEFINITION

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

When adapting the module to another lattice, edit section 1.1.1 only. The calculation and update code below that section should not require lattice-specific changes.

Preparing the Full Ring

Import the module and prepare the lattice with:

import linear_lattice as lin

magnets, lattice, data, correction, parameters = lin.prepare_lattice()

prepare_lattice() performs the complete initialization:

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

The current ring contains:

20 repeated cells
2,220 ordered elements

There is no cell/ring selection. prepare_lattice() always returns the full ring and all ring-level calculations use a repetition factor of one.

Returned Values
Value	Meaning
magnets	Unique magnet-family objects
lattice	Complete ordered ring
data	Linear-optics results and sampled functions
correction	Chromatic correction result, or None if disabled
parameters	Parameter dictionary synchronized with the prepared ring

The chromatic correction has the form:

[
    family1,
    strength1,
    family2,
    strength2,
    corrected_chrom_x,
    corrected_chrom_y,
]

Important linear results can be accessed with:

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

CELL_M4 and CELL_M5 remain available as compatibility aliases for LATTICE_M4 and LATTICE_M5.

Updating the Prepared Ring

After initialization, use update_linear() instead of rebuilding the ring during every optimizer evaluation.

The correct order is:

Change values in the returned parameters dictionary.
Record the exact names of the changed parameters.
Call update_linear().
Use the returned lattice, data, and parameters.

For example:

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

The dictionary contains the new values. edited_variables tells the function which values must be copied into the ring and which calculations are no longer valid.

Changing the dictionary without listing the parameter does not update the magnet:

parameters["X1"] = new_X1

# Incorrect: X1 is missing from edited_variables.
lin.update_linear(lattice, data, parameters, [])

Listing a parameter before changing its dictionary value simply writes the old value back into the ring.

Recalculation Rules
Changed parameters	Magnet values	Linear optics	Chromatic correction
energy, LSD, F, X1–X9	Update	Recompute	Recompute
Sextupole parameters	Update	Reuse	Recompute
ko1, ko2, ko3	Update	Reuse	Reuse

The dependency logic is:

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

When variables from multiple categories are passed together, the function performs every calculation required by the most influential category.

Understanding new_correction

For a nonlinear-only update:

new_correction is None

This means that chromatic correction was not recalculated. It does not mean that the corrected sextupole strengths were removed from the ring.

Preserve the previous correction information with:

lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    edited_variables,
)

if new_correction is not None:
    correction = new_correction
Updating an Optimizer Candidate

An optimizer can provide a dictionary containing one candidate:

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

The nonlinear calculation should always receive the returned lattice and data, because linear-variable changes invalidate the previous optics.

Chromatic Correction Families

The default correction families and targets are defined in section 1.1.1:

CHROMATIC_FAMILY1 = "SF1"
CHROMATIC_FAMILY2 = "SD1"
TARGET_CHROM_X = 0.0
TARGET_CHROM_Y = 0.0

When chromatic correction is enabled, the solved correction-family strengths are written into the ring and synchronized with the parameter dictionary through CORRECTION_PARAMETER_MAP.

If a correction-family strength is supplied as an edited variable, it may immediately be replaced by the strength required to reach the target chromaticities. To preserve an explicitly supplied value, disable the correction for that update:

lattice, data, new_correction, parameters = lin.update_linear(
    lattice,
    data,
    parameters,
    ["ksf1", "ksd1"],
    correct_chromatic=False,
)
Complete Rebuild

Use prepare_lattice() again instead of update_linear() when:

starting an independent calculation,
changing the ring ordering,
changing the number of repeated cells,
adding or removing magnets,
changing magnet-family names,
the parameter dictionary and ring may no longer be synchronized.

For example:

p = lin.PARAMETERS.copy()
p["X1"] = new_X1

magnets, lattice, data, correction, parameters = lin.prepare_lattice(
    parameters=p
)

update_linear() changes supported parameter values. It does not change the ring topology.

Interface With nonlinear.py

The intended initialization is:

import linear_lattice as lin

magnets, lattice, data, correction, parameters = lin.prepare_lattice()

The nonlinear module should use the returned lattice and data rather than reconstructing the accelerator internally.

During optimization:

optimizer candidate
        ↓
update parameters
        ↓
update_linear()
        ↓
updated full ring + valid linear data
        ↓
nonlinear calculation

This keeps the linear and nonlinear models synchronized while avoiding unnecessary optics calculations.

Optional Validation and Plotting

The final section of linear_lattice.py contains optional utilities for:

consistency checks,
terminal summaries,
plotting Twiss functions and dispersion,
standalone execution.

Run the module directly with:

python linear_lattice.py

Run the detailed report with:

python test_linear_lattice.py

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

In short:

prepare once
    → edit the parameter dictionary
    → call update_linear()
    → use the returned ring and linear data

