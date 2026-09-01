## Questions and Further Work

* Review whether `prepare_lattice()` and `update_linear()` return all the necessary information in the most convenient format.
* Verify that lattice updates are correctly propagated to the nonlinear calculations.
* The linear module now supports analyzing either one cell or the full ring. However, the nonlinear and optimization modules may still contain cell-dependent assumptions.
* Compare the nonlinear results obtained after changing linear variables with independent OPA simulations.

## Linear Lattice Report

The following test validates the linear lattice calculations, the selection between one cell and the full ring, and the behavior of `update_linear()`:

```python
import linear_lattice as lin
from contextlib import redirect_stdout
from pathlib import Path


REPORT_FILE = "linear_lattice_report.txt"
SCOPE = "cell"


def run_report():
    repetitions = lin.N_CELLS if SCOPE == "cell" else 1

    magnets, lattice, data, correction, parameters = (
        lin.prepare_lattice(scope=SCOPE)
    )

    print("=" * 80)
    print("LINEAR LATTICE TEST REPORT")
    print("=" * 80)

    print("\nGENERAL")
    print(f"Defined magnets       : {len(magnets)}")
    print(f"Analysis scope        : {SCOPE}")
    print(f"Elements analyzed     : {len(lattice)}")
    print(f"Elements in one cell  : {len(lin.CELL_NAMES)}")
    print(f"Number of cells       : {lin.N_CELLS}")
    print(f"Elements in full ring : {len(lin.RING_NAMES)}")
    print(f"Energy [GeV]          : {parameters['energy']}")

    print("\nPERIODIC TWISS AT s = 0")
    bx, ax, gx, by, ay, gy = data[lin.CS0]

    print(f"beta_x  = {bx:.12g}")
    print(f"alpha_x = {ax:.12g}")
    print(f"gamma_x = {gx:.12g}")
    print(f"beta_y  = {by:.12g}")
    print(f"alpha_y = {ay:.12g}")
    print(f"gamma_y = {gy:.12g}")

    print("\nPERIODIC DISPERSION AT s = 0")
    print(f"Dx  = {data[lin.DISP0][0]:.12g}")
    print(f"Dpx = {data[lin.DISP0][1]:.12g}")

    print("\nRING OPTICS")
    print(f"Tune x                  = {data[lin.TUNE_X]:.12g}")
    print(f"Tune y                  = {data[lin.TUNE_Y]:.12g}")
    print(f"Natural chromaticity x  = {data[lin.CHROM_X]:.12g}")
    print(f"Natural chromaticity y  = {data[lin.CHROM_Y]:.12g}")
    print(f"Natural emittance [m]   = {data[lin.EMITTANCE]:.12g}")
    print(f"Circumference [m]       = {data[lin.CIRCUMFERENCE]:.12g}")

    print("\nCHROMATIC CORRECTION")

    if correction is None:
        print("Chromatic correction disabled.")
    else:
        print(f"{correction[0]} = {correction[1]:.12g}")
        print(f"{correction[2]} = {correction[3]:.12g}")
        print(f"Corrected chromaticity x = {correction[4]:.12g}")
        print(f"Corrected chromaticity y = {correction[5]:.12g}")

    print(f"\n{SCOPE.upper()} MATRICES")
    print("\nM4 =")
    print(data[lin.CELL_M4])
    print("\nM5 =")
    print(data[lin.CELL_M5])

    checks = lin.check_linear_lattice(
        lattice,
        data,
        repetitions=repetitions,
    )

    print("\nCONSISTENCY CHECKS")
    print(f"Symplectic error        = {checks[lin.CHECK_SYMPLECTIC]:.6e}")
    print(f"Twiss closure           = {checks[lin.CHECK_TWISS_CLOSURE]:.6e}")
    print(f"Dispersion closure      = {checks[lin.CHECK_DISPERSION_CLOSURE]:.6e}")
    print(f"CS identity x error     = {checks[lin.CHECK_CS_X]:.6e}")
    print(f"CS identity y error     = {checks[lin.CHECK_CS_Y]:.6e}")
    print(f"Total ring bending [deg]= {checks[lin.CHECK_RING_BENDING_DEG]:.12f}")

    print("\nSAMPLING")
    print(f"Number of s samples     = {len(data[lin.S_VALUES])}")
    print(f"Analyzed length [m]     = {data[lin.S_VALUES][-1]:.12g}")

    print("\nMAGNET TABLE")
    print(
        f"{'NAME':<10} {'TYPE':<12} {'LENGTH':>12} "
        f"{'ANGLE':>12} {'K':>14} {'S':>14} {'O':>14}"
    )
    print("-" * 92)

    for elem in lin.unique_magnets(lattice):
        print(
            f"{elem[lin.NAME]:<10} "
            f"{elem[lin.TYPE]:<12} "
            f"{elem[lin.LENGTH]:>12.6g} "
            f"{elem[lin.ANGLE]:>12.6g} "
            f"{elem[lin.K]:>14.6g} "
            f"{elem[lin.S]:>14.6g} "
            f"{elem[lin.O]:>14.6g}"
        )

    old_data = data
    parameters["ko1"] += 1e-6

    lattice, data, update_correction, parameters = lin.update_linear(
        lattice,
        data,
        parameters,
        ["ko1"],
    )

    nonlinear_ok = data is old_data and update_correction is None

    old_data = data
    parameters["ks1"] += 1e-6

    lattice, data, update_correction, parameters = lin.update_linear(
        lattice,
        data,
        parameters,
        ["ks1"],
    )

    chromatic_ok = data is old_data and update_correction is not None

    old_data = data
    parameters["X1"] += 1e-6

    lattice, data, update_correction, parameters = lin.update_linear(
        lattice,
        data,
        parameters,
        ["X1"],
    )

    linear_ok = data is not old_data and update_correction is not None

    print("\nUPDATE_LINEAR TESTS")
    print(f"Nonlinear-only update  : {'PASS' if nonlinear_ok else 'FAIL'}")
    print(f"Chromatic-only update  : {'PASS' if chromatic_ok else 'FAIL'}")
    print(f"Linear-optics update   : {'PASS' if linear_ok else 'FAIL'}")

    if not all((nonlinear_ok, chromatic_ok, linear_ok)):
        raise AssertionError("One or more update_linear tests failed.")

    print("\n" + "=" * 80)
    print("TEST COMPLETED SUCCESSFULLY")
    print("=" * 80)


def main():
    report_path = Path(REPORT_FILE)

    run_report()

    with report_path.open("w", encoding="utf-8") as file:
        with redirect_stdout(file):
            run_report()

    print(f"\nReport saved to: {report_path.resolve()}")


if __name__ == "__main__":
    main()
```

The test was successfully executed for both `scope="cell"` and `scope="ring"`. The calculated linear-optics quantities agree with the independent OPA results shown in the reference screenshots.

The three update cases also passed:

* Nonlinear-only changes update the magnet strength without recalculating the linear optics.
* Sextupole changes reuse the existing linear-optics data and repeat only the chromatic correction.
* Linear-variable changes recalculate both the linear optics and the chromatic correction.

Therefore, the linear lattice implementation and its update logic have been successfully validated at the linear level. The next step is to verify that the updated lattice is interpreted correctly by the nonlinear module.


