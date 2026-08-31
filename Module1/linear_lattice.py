
######################################33
#MODULE 1: LINEAR OPTICS
"""Linear lattice contains the manually input ESRF-style cell used in this project.

The module contains only the lattice construction, linear transfer matrices,
periodic linear optics, first-order chromaticity, and chromatic correction.

Magnet representation (Python list):
    [name, type, length, angle, K, S, O, M, M5]

Recommended use from a nonlinear program:
    import linear_lattice as lin
    magnets, cell, data, correction, parameters = lin.prepare_lattice()
"""

import math
import numpy as np

# CHAPTER 1.1: USER INPUT (OPA TRANSLATION)

PARAMETERS = {
    "energy": 3.0,
    "LSD": 0.1,   
    "F": 0.8,

    # Variables that may later be changed by an optimizer.
    "X1": 3.633167514008421,
    "X2": -4.258277621861492,
    "X3": -2.690860661253351,
    "X4": 2.754505457254375,
    "X5": -3.336431720192718,
    "X6": -1.492176552197721,
    "X7": 2.943728718874649e-3,
    "X8": 6.010287762639232e-1,
    "X9": 5.370545478298007e-1,

    "kse1": -31.67808026513335,
    "kfd2": -51.10107383195414,
    "kfd3": -93.324704323109,
    "ks1": -31.98668968024124,
    "ks2": 41.64059482880812,
    "ksd3": -11.2882623813059,
    "ks1s": 696.35209098583,
    "ks2s": 488.0292055866639,
    #I need to remove the length of 10^-8
    "ko1": 4.72312991538625,
    "ko2": 46.6898786404116,
    "ko3": 23.51531647572548,
    "ksf1": 29.726909817,
    "ksd1": -127.01,
}

VARY = ["X1", "X2", "X3", "X4", "X5", "X6"]



# 1.1.2 MAGNET LIST FORMAT
# Every magnet is: [name, type, length, angle, K, S, O, M, M5]

NAME = 0;TYPE = 1; LENGTH = 2; ANGLE = 3; 
K = 4; S = 5; O = 6; M = 7; M5 = 8


def magnet(name, magnet_type, length, angle=0.0, K_value=0.0, S_value=0.0, O_value=0.0):
    return [
        name,
        magnet_type,
        float(length),
        float(angle),
        float(K_value),
        float(S_value),
        float(O_value),
        None, #Free if needed?
        None,
    ]



# 1.1.3. LINEAR TRANSFER MATRICES


def plane_matrix(L, k):
    """Generic 2x2 uncoupled linear map for constant focusing k."""
    if abs(k) < 1e-15:
        return np.array([[1.0, L], [0.0, 1.0]])

    if k > 0.0:
        root = math.sqrt(k)
        C = math.cos(root * L)
        Sine = math.sin(root * L) / root
    else:
        root = math.sqrt(-k)
        C = math.cosh(root * L)
        Sine = math.sinh(root * L) / root

    return np.array([[C, Sine], [-k * Sine, C]])


def drift4(L):
    matrix = np.eye(4)
    matrix[0:2, 0:2] = plane_matrix(L, 0.0)
    matrix[2:4, 2:4] = plane_matrix(L, 0.0)
    return matrix


def drift5(L):
    matrix = np.eye(5)
    matrix[:4, :4] = drift4(L)
    return matrix


def quadrupole4(L, k):
    matrix = np.zeros((4, 4))
    matrix[0:2, 0:2] = plane_matrix(L, k)
    matrix[2:4, 2:4] = plane_matrix(L, -k)
    return matrix


def quadrupole5(L, k):
    matrix = np.eye(5)
    matrix[:4, :4] = quadrupole4(L, k)
    return matrix


def bending4(L, k, h):
    matrix = np.zeros((4, 4))
    matrix[0:2, 0:2] = plane_matrix(L, k + h**2)
    matrix[2:4, 2:4] = plane_matrix(L, -k)
    return matrix


def bending5(L, k, h):
    matrix = np.eye(5)
    matrix[:4, :4] = bending4(L, k, h)

    kx = k + h**2
    Cx = matrix[0, 0]
    Sx = matrix[0, 1]

    if abs(kx) < 1e-15:
        matrix[0, 4] = 0.5 * h * L**2
    else:
        matrix[0, 4] = h * (1.0 - Cx) / kx

    matrix[1, 4] = h * Sx
    return matrix


def transfer4(elem, L=None):
    """4x4 map through a complete magnet or through a piece of it."""
    ds = elem[LENGTH] if L is None else float(L)

    if elem[TYPE] in ("drift", "sextupole"):
        return drift4(ds)

    if elem[TYPE] == "quadrupole":
        return quadrupole4(ds, elem[K])

    if elem[TYPE] == "bending":
        if elem[LENGTH] == 0.0:
            return np.eye(4)
        h = math.radians(elem[ANGLE]) / elem[LENGTH]
        return bending4(ds, elem[K], h)

    if elem[TYPE] == "multipole":
        return np.eye(4)

    raise ValueError("Unknown magnet type: " + str(elem[TYPE]))


def transfer5(elem, L=None):
    """5x5 map used to propagate horizontal dispersion."""
    ds = elem[LENGTH] if L is None else float(L)

    if elem[TYPE] in ("drift", "sextupole"):
        return drift5(ds)

    if elem[TYPE] == "quadrupole":
        return quadrupole5(ds, elem[K])

    if elem[TYPE] == "bending":
        if elem[LENGTH] == 0.0:
            return np.eye(5)
        h = math.radians(elem[ANGLE]) / elem[LENGTH]
        return bending5(ds, elem[K], h)

    if elem[TYPE] == "multipole":
        return np.eye(5)

    raise ValueError("Unknown magnet type: " + str(elem[TYPE]))


def compute_matrices(elem):
    """Compute and store the full linear maps of one magnet."""
    elem[M] = transfer4(elem)
    elem[M5] = transfer5(elem)


def unique_magnets(lattice):
    """Return each magnet object in a lattice only once, preserving order."""
    result = []
    seen = set()
    for elem in lattice:
        marker = id(elem)
        if marker not in seen:
            seen.add(marker)
            result.append(elem)
    return result


def refresh_matrices(lattice):
    """Recompute M and M5 for all unique magnets in a lattice."""
    for elem in unique_magnets(lattice):
        compute_matrices(elem)


# 1.1.4. MAGNET DEFINITIONS


def define_magnets(p):
    LSD = p["LSD"]
    F = p["F"]

    magnets = [
        # Drifts
        magnet("D1", "drift", 2.654400 - LSD),
        magnet("D4", "drift", 0.081240),
        magnet("D11", "drift", 0.063628),
        magnet("D12", "drift", 0.0099526),
        magnet("D5D6", "drift", p["X8"]),
        magnet("D9D10", "drift", p["X9"]),

        # Quadrupoles
        magnet("QF1", "quadrupole", 0.349140, K_value=p["X1"]),
        magnet("QD2", "quadrupole", 0.222950, K_value=p["X2"]),
        magnet("QD3", "quadrupole", 0.194780, K_value=p["X3"]),
        magnet("QF4", "quadrupole", 0.224580, K_value=p["X4"]),
        magnet("QD5", "quadrupole", 0.210950, K_value=p["X5"]),
        magnet("QF7", "quadrupole", 0.020986, K_value=p["X6"]),

        # Sextupoles: on-momentum linear map = drift
        magnet("SE1", "sextupole", LSD, S_value=p["kse1"]),
        magnet("FD2", "sextupole", 0.094502, S_value=p["kfd2"]),
        magnet("FD3", "sextupole", p["X7"], S_value=p["kfd3"]),
        magnet("S1", "sextupole", LSD, S_value=p["ks1"]),
        magnet("S2", "sextupole", LSD, S_value=p["ks2"]),
        magnet("SD3", "sextupole", 0.010176, S_value=p["ksd3"]),
        magnet("S1S", "sextupole", 0.002964, S_value=p["ks1s"]),
        magnet("S2S", "sextupole", 0.172130, S_value=p["ks2s"]),
        magnet("SF1", "sextupole", 0.220440, S_value=p["ksf1"]),
        magnet("SD1", "sextupole", LSD, S_value=p["ksd1"]),

        # Thin higher multipoles: O stores the integrated nonlinear strength.
        # They have zero physical length and therefore no linear effect.
        magnet("O1", "multipole", 0.0, O_value=p["ko1"]),  #cero values are going to be treated differently
        magnet("O2", "multipole", 0.0, O_value=p["ko2"]),
        magnet("O3", "multipole", 0.0, O_value=p["ko3"]),

        # Bending / combined-function magnets
        magnet("DQ6", "bending", 0.275390, angle=-0.73179259 * F, K_value=2.692600),
        magnet("A1", "bending", 0.075497, angle=0.0021719 * F),
        magnet("A2", "bending", 0.384040, angle=0.53380 * F),
        magnet("A3", "bending", 0.001995, angle=0.00032534 * F),
        magnet("A4", "bending", 0.913400, angle=2.0382 * F),
        magnet("A5", "bending", 0.152490, angle=0.93133 * F),
        magnet("B1", "bending", 0.400570, angle=0.63294 * F),
        magnet("B2", "bending", 0.563170, angle=1.1254 * F),
        magnet("B3", "bending", 0.362720, angle=1.1741 * F),
        magnet("B4", "bending", 0.285610, angle=1.4465 * F),
        magnet("B5", "bending", 0.240960, angle=0.58358 * F),
        magnet("B1S", "bending", 0.015767, angle=0.080780 * F),
        magnet("B2S", "bending", 0.001644, angle=-0.00041155 * F),
        magnet("B3S", "bending", 0.212550, angle=1.7586 * F),
        magnet("DQ1S", "bending", 0.257080, angle=0.81690 * F, K_value=-5.135300),
        magnet("ABQ1", "bending", 0.215990, angle=-0.60542 * F, K_value=6.191000),
    ]

    for elem in magnets:
        compute_matrices(elem)

    return magnets



# 1.1.5 MANUAL LATTICE


DA1 = ["A1", "A2", "A3", "A4", "A5"]
IDA1 = DA1[::-1]

DBA = [
    "D1", "SE1", "QF1", "FD2", "QD2", "FD3",
    *IDA1,
    "D4", "QD3", "SD1", "O2", "D5D6", "S1", "QF4", "SF1",
    "O1", "QF4", "S2", "D9D10", "O3", "SD1", "QD5", "D11",
    "B1", "B2", "B3", "B4", "B5", "D12", "QF7", "SD3", "DQ6",
]

CELA = [
    "S1S", "ABQ1", "S2S", "DQ1S", "B1S", "B2S", "B3S",
    "B2S", "B1S", "DQ1S", "S2S", "ABQ1", "S1S",
]

CELL_NAMES = DBA + CELA + CELA + CELA + DBA[::-1]
N_CELLS = 20
RING_NAMES = CELL_NAMES * N_CELLS


def build_lattice(names, magnets):
    """Replace every lattice name by its magnet list."""
    magnet_by_name = {elem[NAME]: elem for elem in magnets}

    missing = [name for name in names if name not in magnet_by_name]
    if missing:
        raise KeyError("Undefined magnets: " + str(sorted(set(missing))))

    return [magnet_by_name[name] for name in names]


def get_magnet(lattice, magnet_name):
    """Return the first magnet with the requested name."""
    for elem in lattice:
        if elem[NAME] == magnet_name:
            return elem
    raise KeyError(f"Magnet {magnet_name} was not found in the lattice.")


def set_quadrupole_strength(lattice, magnet_name, strength):
    """Update K and immediately refresh the corresponding linear matrices."""
    elem = get_magnet(lattice, magnet_name)
    if elem[TYPE] != "quadrupole":
        raise ValueError(f"{magnet_name} is not a quadrupole.")
    elem[K] = float(strength)
    compute_matrices(elem)


def set_sextupole_strength(lattice, family_name, strength):
    """Update S. No M/M5 refresh is needed because its linear map is a drift."""
    elem = get_magnet(lattice, family_name)
    if elem[TYPE] != "sextupole":
        raise ValueError(f"{family_name} is not a sextupole.")
    elem[S] = float(strength)

# 1.1.6 COURANT-SNYDER TRANSPORT

def courant_snyder_matrix(matrix):
    return np.array([
        [matrix[0, 0]**2, -2*matrix[0, 0]*matrix[0, 1], matrix[0, 1]**2, 0, 0, 0],
        [-matrix[0, 0]*matrix[1, 0], matrix[0, 0]*matrix[1, 1] + matrix[0, 1]*matrix[1, 0], -matrix[0, 1]*matrix[1, 1], 0, 0, 0],
        [matrix[1, 0]**2, -2*matrix[1, 0]*matrix[1, 1], matrix[1, 1]**2, 0, 0, 0],
        [0, 0, 0, matrix[2, 2]**2, -2*matrix[2, 2]*matrix[2, 3], matrix[2, 3]**2],
        [0, 0, 0, -matrix[2, 2]*matrix[3, 2], matrix[2, 2]*matrix[3, 3] + matrix[2, 3]*matrix[3, 2], -matrix[2, 3]*matrix[3, 3]],
        [0, 0, 0, matrix[3, 2]**2, -2*matrix[3, 2]*matrix[3, 3], matrix[3, 3]**2],
    ])


def propagate_linear_functions(lattice, cs0, disp0, step=0.01):
    s = 0.0
    mux = 0.0
    muy = 0.0
    chromx = 0.0
    chromy = 0.0
    radiation = np.zeros(6)

    cs = cs0.copy()
    disp = disp0.copy()

    s_values = [s]
    cs_values = [cs.copy()]
    disp_values = [disp.copy()]

    for elem in lattice:
        if elem[LENGTH] == 0.0:
            continue

        n_full = int(elem[LENGTH] // step)
        pieces = [step] * n_full
        remainder = elem[LENGTH] - n_full * step

        if remainder > 1e-14:
            pieces.append(remainder)

        for ds in pieces:
            beta_x_i = cs[0]
            beta_y_i = cs[3]

            if elem[K] != 0.0:
                chromx += ds * beta_x_i * elem[K]
                chromy += ds * beta_y_i * elem[K]

            if elem[TYPE] == "bending" and elem[ANGLE] != 0.0:
                h = math.radians(elem[ANGLE]) / elem[LENGTH]
                Hx = (
                    beta_x_i * disp[1]**2
                    + 2.0 * cs[1] * disp[0] * disp[1]
                    + cs[2] * disp[0]**2
                )

                radiation[0] += ds * disp[0] * h
                radiation[1] += ds * h**2
                radiation[2] += ds * abs(h)**3
                radiation[3] += ds * (1.0 + 2.0 * elem[K] / h**2) * disp[0] * h**3
                radiation[4] += ds * Hx * abs(h)**3
                radiation[5] += ds * (disp[0] * elem[K])**2

            matrix4 = transfer4(elem, ds)
            matrix5 = transfer5(elem, ds)

            cs = courant_snyder_matrix(matrix4) @ cs
            disp = matrix5 @ disp

            beta_x_f = cs[0]
            beta_y_f = cs[3]

            argx = matrix4[0, 1] / math.sqrt(beta_x_i * beta_x_f)
            argy = matrix4[2, 3] / math.sqrt(beta_y_i * beta_y_f)

            mux += math.asin(np.clip(argx, -1.0, 1.0))
            muy += math.asin(np.clip(argy, -1.0, 1.0))

            s += ds
            s_values.append(s)
            cs_values.append(cs.copy())
            disp_values.append(disp.copy())

    return [
        mux,
        muy,
        np.array(s_values),
        np.array(cs_values),
        np.array(disp_values),
        radiation,
        chromx,
        chromy,
    ]




def lattice_matrix(lattice):
    matrix4 = np.eye(4)
    matrix5 = np.eye(5)

    for elem in lattice:
        matrix4 = elem[M] @ matrix4
        matrix5 = elem[M5] @ matrix5

    return matrix4, matrix5


def periodic_twiss_and_dispersion(lattice):
    matrix4, matrix5 = lattice_matrix(lattice)

    trace_x = matrix4[0, 0] + matrix4[1, 1]
    trace_y = matrix4[2, 2] + matrix4[3, 3]

    if abs(trace_x) >= 2.0 or abs(trace_y) >= 2.0:
        raise ValueError(
            f"Unstable lattice: horizontal trace={trace_x:.8f}, vertical trace={trace_y:.8f}"
        )

    sin_mux = np.sign(matrix4[0, 1]) * math.sqrt(
        -matrix4[0, 1] * matrix4[1, 0]
        - (matrix4[0, 0] - matrix4[1, 1])**2 / 4.0
    )

    sin_muy = np.sign(matrix4[2, 3]) * math.sqrt(
        -matrix4[2, 3] * matrix4[3, 2]
        - (matrix4[2, 2] - matrix4[3, 3])**2 / 4.0
    )

    ax = (matrix4[0, 0] - matrix4[1, 1]) / (2.0 * sin_mux)
    bx = matrix4[0, 1] / sin_mux
    gx = (1.0 + ax**2) / bx

    ay = (matrix4[2, 2] - matrix4[3, 3]) / (2.0 * sin_muy)
    by = matrix4[2, 3] / sin_muy
    gy = (1.0 + ay**2) / by

    denominator = 2.0 - matrix5[0, 0] - matrix5[1, 1]

    disp = (
        matrix5[0, 1] * matrix5[1, 4]
        + matrix5[0, 4] * (1.0 - matrix5[1, 1])
    ) / denominator

    dispd = (
        matrix5[1, 0] * matrix5[0, 4]
        + matrix5[1, 4] * (1.0 - matrix5[0, 0])
    ) / denominator

    cs0 = np.array([bx, ax, gx, by, ay, gy])
    disp0 = np.array([disp, dispd, 0.0, 0.0, 1.0])

    return cs0, disp0, matrix4, matrix5


def linear_optics(cell, energy, repetitions=20, step=0.01):
    cs0, disp0, matrix4, matrix5 = periodic_twiss_and_dispersion(cell)

    mux, muy, s_values, cs_values, disp_values, radiation, chromx_i, chromy_i = (
        propagate_linear_functions(cell, cs0, disp0, step)
    )

    I2 = radiation[1]
    I4 = radiation[3]
    I5 = radiation[4]

    if I2 == 0.0 or I2 - I4 == 0.0:
        natural_emittance = float("nan")
    else:
        natural_emittance = (
            3.8319e-13
            * (1000.0 * energy / 0.5109989)**2
            * I5 / (I2 - I4)
        )

    tune_x = repetitions * mux / (2.0 * math.pi)
    tune_y = repetitions * muy / (2.0 * math.pi)
    chrom_x = -(repetitions / (4.0 * math.pi)) * chromx_i
    chrom_y = +(repetitions / (4.0 * math.pi)) * chromy_i
    circumference = repetitions * s_values[-1]

    return [
        cs0,
        disp0,
        matrix4,
        matrix5,
        tune_x,
        tune_y,
        chrom_x,
        chrom_y,
        natural_emittance,
        circumference,
        s_values,
        cs_values,
        disp_values,
    ]


# Indices for the linear_optics result list.
CS0 = 0
DISP0 = 1
CELL_M4 = 2
CELL_M5 = 3
TUNE_X = 4
TUNE_Y = 5
CHROM_X = 6
CHROM_Y = 7
EMITTANCE = 8
CIRCUMFERENCE = 9
S_VALUES = 10
CS_VALUES = 11
DISP_VALUES = 12



# 1.1.8. FIRST-ORDER CHROMATIC CORRECTION


def chromaticity_with_sextupoles(lattice, repetitions=20, step=0.01):
    """First-order chromaticities using K_eff = K - 2*Dx*S."""
    cs, disp, _, _ = periodic_twiss_and_dispersion(lattice)

    integral_x = 0.0
    integral_y = 0.0

    for elem in lattice:
        if elem[LENGTH] == 0.0:
            continue

        n_full = int(elem[LENGTH] // step)
        pieces = [step] * n_full
        remainder = elem[LENGTH] - n_full * step

        if remainder > 1e-14:
            pieces.append(remainder)

        for ds in pieces:
            k_eff = elem[K] - 2.0 * disp[0] * elem[S]

            integral_x += ds * cs[0] * k_eff
            integral_y += ds * cs[3] * k_eff

            cs = courant_snyder_matrix(transfer4(elem, ds)) @ cs
            disp = transfer5(elem, ds) @ disp

    chrom_x = -(repetitions / (4.0 * math.pi)) * integral_x
    chrom_y = +(repetitions / (4.0 * math.pi)) * integral_y

    return chrom_x, chrom_y


def correct_chromaticity(
    lattice,
    family1="SF1",
    family2="SD1",
    target_x=0.0,
    target_y=0.0,
    repetitions=20,
    step=0.01,
):
    """Solve the two-family first-order chromatic correction problem."""
    if family1 == family2:
        raise ValueError("family1 and family2 must be different sextupole families.")

    cs, disp, _, _ = periodic_twiss_and_dispersion(lattice)

    base_x = 0.0
    base_y = 0.0
    response_x1 = 0.0
    response_y1 = 0.0
    response_x2 = 0.0
    response_y2 = 0.0

    found1 = False
    found2 = False

    for elem in lattice:
        if elem[LENGTH] == 0.0:
            continue

        if elem[NAME] == family1:
            if elem[TYPE] != "sextupole":
                raise ValueError(f"{family1} is not a sextupole.")
            found1 = True

        if elem[NAME] == family2:
            if elem[TYPE] != "sextupole":
                raise ValueError(f"{family2} is not a sextupole.")
            found2 = True

        n_full = int(elem[LENGTH] // step)
        pieces = [step] * n_full
        remainder = elem[LENGTH] - n_full * step

        if remainder > 1e-14:
            pieces.append(remainder)

        for ds in pieces:
            beta_x = cs[0]
            beta_y = cs[3]
            Dx = disp[0]

            if elem[NAME] == family1:
                response_x1 += 2.0 * ds * beta_x * Dx
                response_y1 += 2.0 * ds * beta_y * Dx
                S_here = 0.0
            elif elem[NAME] == family2:
                response_x2 += 2.0 * ds * beta_x * Dx
                response_y2 += 2.0 * ds * beta_y * Dx
                S_here = 0.0
            else:
                S_here = elem[S]

            k_eff = elem[K] - 2.0 * Dx * S_here
            base_x += ds * beta_x * k_eff
            base_y += ds * beta_y * k_eff

            cs = courant_snyder_matrix(transfer4(elem, ds)) @ cs
            disp = transfer5(elem, ds) @ disp

    if not found1:
        raise KeyError(f"Sextupole family {family1} was not found in the lattice.")
    if not found2:
        raise KeyError(f"Sextupole family {family2} was not found in the lattice.")

    A = np.array([
        [response_x1, response_x2],
        [response_y1, response_y2],
    ])

    b = np.array([
        base_x + target_x * (4.0 * math.pi / repetitions),
        base_y - target_y * (4.0 * math.pi / repetitions),
    ])

    if abs(np.linalg.det(A)) < 1e-14:
        raise ValueError(
            "The two sextupole families do not provide independent chromatic correction."
        )

    S1, S2 = np.linalg.solve(A, b)

    set_sextupole_strength(lattice, family1, S1)
    set_sextupole_strength(lattice, family2, S2)

    corrected_x, corrected_y = chromaticity_with_sextupoles(
        lattice,
        repetitions=repetitions,
        step=step,
    )

    return S1, S2, corrected_x, corrected_y


def correct_chromaticity_from_data(
    lattice,
    data,
    family1="SF1",
    family2="SD1",
    target_x=0.0,
    target_y=0.0,
    repetitions=20,
    step=0.01,
):
    """Chromatic correction reusing optics already computed by linear_optics().

    This avoids recomputing periodic Twiss/dispersion and avoids a second full
    linear propagation during every optimizer evaluation.  The sampled
    cs/disp arrays in ``data`` are exactly the values needed before each step.
    """
    if family1 == family2:
        raise ValueError("family1 and family2 must be different sextupole families.")

    cs_values = np.asarray(data[CS_VALUES], dtype=float)
    disp_values = np.asarray(data[DISP_VALUES], dtype=float)

    base_x = base_y = 0.0
    response_x1 = response_y1 = 0.0
    response_x2 = response_y2 = 0.0
    found1 = found2 = False
    sample = 0

    for elem in lattice:
        if elem[LENGTH] == 0.0:
            continue

        if elem[NAME] == family1:
            if elem[TYPE] != "sextupole":
                raise ValueError(f"{family1} is not a sextupole.")
            found1 = True
        if elem[NAME] == family2:
            if elem[TYPE] != "sextupole":
                raise ValueError(f"{family2} is not a sextupole.")
            found2 = True

        n_full = int(elem[LENGTH] // step)
        pieces = [step] * n_full
        remainder = elem[LENGTH] - n_full * step
        if remainder > 1e-14:
            pieces.append(remainder)

        for ds in pieces:
            cs = cs_values[sample]
            disp = disp_values[sample]
            beta_x = cs[0]
            beta_y = cs[3]
            Dx = disp[0]

            if elem[NAME] == family1:
                response_x1 += 2.0 * ds * beta_x * Dx
                response_y1 += 2.0 * ds * beta_y * Dx
                S_here = 0.0
            elif elem[NAME] == family2:
                response_x2 += 2.0 * ds * beta_x * Dx
                response_y2 += 2.0 * ds * beta_y * Dx
                S_here = 0.0
            else:
                S_here = elem[S]

            k_eff = elem[K] - 2.0 * Dx * S_here
            base_x += ds * beta_x * k_eff
            base_y += ds * beta_y * k_eff
            sample += 1

    if sample != len(cs_values) - 1:
        raise RuntimeError(
            f"Linear sampling mismatch: used {sample} steps, "
            f"but data contains {len(cs_values) - 1}."
        )
    if not found1:
        raise KeyError(f"Sextupole family {family1} was not found in the lattice.")
    if not found2:
        raise KeyError(f"Sextupole family {family2} was not found in the lattice.")

    A = np.array([
        [response_x1, response_x2],
        [response_y1, response_y2],
    ])
    b = np.array([
        base_x + target_x * (4.0 * math.pi / repetitions),
        base_y - target_y * (4.0 * math.pi / repetitions),
    ])

    if abs(np.linalg.det(A)) < 1e-14:
        raise ValueError(
            "The two sextupole families do not provide independent chromatic correction."
        )

    S1, S2 = np.linalg.solve(A, b)
    set_sextupole_strength(lattice, family1, S1)
    set_sextupole_strength(lattice, family2, S2)

    # No extra lattice traversal is necessary.  Evaluate the corrected
    # first-order chromaticities from the same integrals used in the solve.
    corrected_integral_x = base_x - response_x1 * S1 - response_x2 * S2
    corrected_integral_y = base_y - response_y1 * S1 - response_y2 * S2
    corrected_x = -(repetitions / (4.0 * math.pi)) * corrected_integral_x
    corrected_y = +(repetitions / (4.0 * math.pi)) * corrected_integral_y

    return S1, S2, corrected_x, corrected_y



# 1.1.9. PUBLIC PREPARATION FUNCTION FOR THE NONLINEAR MODULE


def prepare_lattice(
    parameters=None,
    correct_chromatic=True,
    family1="SF1",
    family2="SD1",
    target_chrom_x=0.0,
    target_chrom_y=0.0,
    repetitions=N_CELLS,
    step=0.01,
):
    """Build and prepare the lattice for later nonlinear calculations.

    """
    p = PARAMETERS.copy() if parameters is None else parameters.copy()

    magnets = define_magnets(p)
    cell = build_lattice(CELL_NAMES, magnets)

    data = linear_optics(
        cell,
        energy=p["energy"],
        repetitions=repetitions,
        step=step,
    )

    correction = None

    if correct_chromatic:
        S1, S2, corrected_x, corrected_y = correct_chromaticity_from_data(
            cell,
            data,
            family1=family1,
            family2=family2,
            target_x=target_chrom_x,
            target_y=target_chrom_y,
            repetitions=repetitions,
            step=step,
        )

        correction = [family1, S1, family2, S2, corrected_x, corrected_y]

        # Keep the returned parameter copy synchronized for the default families.
        if family1 == "SF1":
            p["ksf1"] = S1
        elif family1 == "SD1":
            p["ksd1"] = S1

        if family2 == "SF1":
            p["ksf1"] = S2
        elif family2 == "SD1":
            p["ksd1"] = S2

    return magnets, cell, data, correction, p


# 1.1.10. CHECKS / DEBUGGING


def check_linear_lattice(cell, data, repetitions=N_CELLS):
    """Return basic consistency errors for the current linear lattice."""
    M4 = data[CELL_M4]

    J2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
    J4 = np.zeros((4, 4))
    J4[:2, :2] = J2
    J4[2:, 2:] = J2

    symplectic_error = np.linalg.norm(M4.T @ J4 @ M4 - J4)

    cs_values = data[CS_VALUES]
    disp_values = data[DISP_VALUES]

    twiss_closure = np.max(np.abs(cs_values[-1] - cs_values[0]))
    dispersion_closure = np.max(np.abs(disp_values[-1] - disp_values[0]))

    cs_identity_x = np.max(np.abs(cs_values[:, 0] * cs_values[:, 2] - cs_values[:, 1]**2 - 1.0))
    cs_identity_y = np.max(np.abs(cs_values[:, 3] * cs_values[:, 5] - cs_values[:, 4]**2 - 1.0))

    bending_cell = sum(elem[ANGLE] for elem in cell if elem[TYPE] == "bending")
    bending_ring = repetitions * bending_cell

    return [
        symplectic_error,
        twiss_closure,
        dispersion_closure,
        cs_identity_x,
        cs_identity_y,
        bending_ring,
    ]


CHECK_SYMPLECTIC = 0
CHECK_TWISS_CLOSURE = 1
CHECK_DISPERSION_CLOSURE = 2
CHECK_CS_X = 3
CHECK_CS_Y = 4
CHECK_RING_BENDING_DEG = 5


# 1.1.11. OPTIONAL DISPLAY HELPERS


def print_linear_summary(data, energy, correction=None):
    bx, ax, gx, by, ay, gy = data[CS0]

    print("=" * 80)
    print("Beta functions at s = 0")
    print(f"Ax = {ax:.10g}   Ay = {ay:.10g}")
    print(f"Bx = {bx:.10g}   By = {by:.10g}")
    print(f"Gx = {gx:.10g}   Gy = {gy:.10g}")
    print("=" * 80)
    print("Ring")
    print(f"Energy          = {energy}")
    print(f"Nux             = {data[TUNE_X]}")
    print(f"Nuy             = {data[TUNE_Y]}")
    print(f"Natural Chromx  = {data[CHROM_X]}")
    print(f"Natural Chromy  = {data[CHROM_Y]}")
    print(f"Emitx           = {data[EMITTANCE]}")
    print(f"Circumference   = {data[CIRCUMFERENCE]}")

    if correction is not None:
        print("-" * 80)
        print("Chromatic correction")
        print(f"{correction[0]} = {correction[1]}")
        print(f"{correction[2]} = {correction[3]}")
        print(f"Corrected Chromx = {correction[4]}")
        print(f"Corrected Chromy = {correction[5]}")

    print("=" * 80)


def plot_linear_functions(data):
    # Imported only when plotting is actually requested.
    import matplotlib.pyplot as plt

    s_values = data[S_VALUES]
    cs_values = data[CS_VALUES]
    disp_values = data[DISP_VALUES]

    plt.plot(s_values, cs_values[:, 0], label=r"$\beta_x$")
    plt.plot(s_values, cs_values[:, 3], label=r"$\beta_y$")
    plt.plot(s_values, 100.0 * disp_values[:, 0], label=r"$100D_x$")
    plt.xlabel("s [m]")
    plt.ylabel("Linear functions [m]")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# Public names intended to be useful to the later nonlinear module.
__all__ = [
    # parameters and lattice definition
    "PARAMETERS", "VARY", "CELL_NAMES", "RING_NAMES", "N_CELLS",
    # magnet indices
    "NAME", "TYPE", "LENGTH", "ANGLE", "K", "S", "O", "M", "M5",
    # construction / updates
    "magnet", "define_magnets", "build_lattice", "get_magnet",
    "unique_magnets", "refresh_matrices", "set_quadrupole_strength",
    "set_sextupole_strength",
    # transfer maps
    "plane_matrix", "drift4", "drift5", "quadrupole4", "quadrupole5",
    "bending4", "bending5", "transfer4", "transfer5", "compute_matrices",
    # linear optics
    "courant_snyder_matrix", "propagate_linear_functions", "lattice_matrix",
    "periodic_twiss_and_dispersion", "linear_optics",
    # linear output indices
    "CS0", "DISP0", "CELL_M4", "CELL_M5", "TUNE_X", "TUNE_Y",
    "CHROM_X", "CHROM_Y", "EMITTANCE", "CIRCUMFERENCE",
    "S_VALUES", "CS_VALUES", "DISP_VALUES",
    # chromatic correction
    "chromaticity_with_sextupoles", "correct_chromaticity",
    "correct_chromaticity_from_data",
    # ready-to-use frontend
    "prepare_lattice",
    # checks
    "check_linear_lattice", "CHECK_SYMPLECTIC", "CHECK_TWISS_CLOSURE",
    "CHECK_DISPERSION_CLOSURE", "CHECK_CS_X", "CHECK_CS_Y",
    "CHECK_RING_BENDING_DEG",
    # optional display
    "print_linear_summary", "plot_linear_functions",
]


def main():
    """Standalone check. This is NOT executed when the module is imported."""
    magnets, cell, data, correction, parameters = prepare_lattice()

    print(f"Defined magnets : {len(magnets)}")
    print(f"Elements in CELL: {len(CELL_NAMES)}")
    print(f"Elements in RING: {len(RING_NAMES)}")
    print_linear_summary(data, parameters["energy"], correction)

    checks = check_linear_lattice(cell, data)
    print("Linear consistency checks")
    print(f"Symplectic error       = {checks[CHECK_SYMPLECTIC]:.6e}")
    print(f"Twiss closure          = {checks[CHECK_TWISS_CLOSURE]:.6e}")
    print(f"Dispersion closure     = {checks[CHECK_DISPERSION_CLOSURE]:.6e}")
    print(f"CS identity x error    = {checks[CHECK_CS_X]:.6e}")
    print(f"CS identity y error    = {checks[CHECK_CS_Y]:.6e}")
    print(f"Total bending [deg]    = {checks[CHECK_RING_BENDING_DEG]:.12f}")


if __name__ == "__main__":
    main()
