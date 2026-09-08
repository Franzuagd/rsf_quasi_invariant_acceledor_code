import subprocess
import re
import os

import multiprocessing
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import time
import scipy
import random
import pygad
from typing import List
from functools import partial
from scipy.stats import skew
from numpy import linalg as LA
from scipy import linalg
from scipy.interpolate import make_interp_spline
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
random.seed(42)
import sympy as sp
import numpy as np
import sympy as sp
from math import comb
import scipy as sc
import scipy.sparse as sps
from collections import defaultdict

np.set_printoptions(
    precision=17,
    suppress=False,
    linewidth=100000,
    threshold=np.inf
)

delta, x, y, px, py = sp.symbols('delta x y px py')
b1, b2, b3, b4, b5 = sp.symbols('b1 b2 b3 b4 b5')
vars = [delta, x, y, px, py]


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
    "ko1": 4.72312991538625,
    "ko2": 46.6898786404116,
    "ko3": 23.51531647572548,
    "ksf1": 29.726909817,
    "ksd1": -127.01,
}

VARY = ["X1", "X2", "X3", "X4", "X5", "X6"]


# ============================================================================
# 2. MAGNET LIST FORMAT
# ============================================================================
# Every magnet is just:
#
# [name, type, length, angle, K, S, O, M, M5]
#
# The constants below avoid mysterious numbers such as magnet[4].

NAME = 0
TYPE = 1
LENGTH = 2
ANGLE = 3
K = 4
S = 5
O = 6
M = 7
M5 = 8


def magnet(name, magnet_type, length, angle=0.0, K_value=0.0, S_value=0.0, O_value=0.0):
    return [name, magnet_type, float(length), float(angle), float(K_value),
            float(S_value), float(O_value), None, None]


# ============================================================================
# 3. LINEAR TRANSFER MATRICES
# ============================================================================

def plane_matrix(L, k):
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

    if elem[TYPE] == "drift" or elem[TYPE] == "sextupole":
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
    """5x5 map including horizontal dispersion."""
    ds = elem[LENGTH] if L is None else float(L)

    if elem[TYPE] == "drift" or elem[TYPE] == "sextupole":
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
    """Store only the full matrices needed later."""
    elem[M] = transfer4(elem)
    elem[M5] = transfer5(elem)



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

        # Sextupoles: linear map = drift
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

        # Thin higher multipoles: no linear effect
        magnet("O1", "multipole", 0.0, O_value=p["ko1"]),
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


# ============================================================================
# 5. MANUAL LATTICE
# ============================================================================
# The lattice itself is now only lists of magnet names.

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
    """Replace each lattice name by its magnet list."""
    magnet_by_name = {elem[NAME]: elem for elem in magnets}

    missing = [name for name in names if name not in magnet_by_name]
    if missing:
        raise KeyError("Undefined magnets: " + str(sorted(set(missing))))

    return [magnet_by_name[name] for name in names]


# ============================================================================
# 6. COURANT-SNYDER TRANSPORT
# ============================================================================

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


# ============================================================================
# 7. PERIODIC LINEAR OPTICS
# ============================================================================

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


def print_linear_summary(data, energy):
    bx, ax, gx, by, ay, gy = data[CS0]

    print("=" * 80)
    print("Beta functions at s = 0")
    print(f"Ax = {ax:.10g}   Ay = {ay:.10g}")
    print(f"Bx = {bx:.10g}   By = {by:.10g}")
    print(f"Gx = {gx:.10g}   Gy = {gy:.10g}")
    print("=" * 80)
    print("Ring")
    print(f"Energy        = {energy}")
    print(f"Nux           = {data[TUNE_X]}")
    print(f"Nuy           = {data[TUNE_Y]}")
    print(f"Chromx        = {data[CHROM_X]}")
    print(f"Chromy        = {data[CHROM_Y]}")
    print(f"Emitx         = {data[EMITTANCE]}")
    print(f"Circumference = {data[CIRCUMFERENCE]}")
    print("=" * 80)


def plot_linear_functions(data):
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


# ============================================================================
# 8. MAIN
# ============================================================================

def main():
    magnets = define_magnets(PARAMETERS)
    cell = build_lattice(CELL_NAMES, magnets)

    print(f"Defined magnets : {len(magnets)}")
    print(f"Elements in CELL: {len(CELL_NAMES)}")
    print(f"Elements in RING: {len(RING_NAMES)}")

    data = linear_optics(
        cell,
        energy=PARAMETERS["energy"],
        repetitions=N_CELLS,
        step=0.01,
    )

    print_linear_summary(data, PARAMETERS["energy"])

    PLOT_LINEAR = True
    if PLOT_LINEAR:
        plot_linear_functions(data)


if __name__ == "__main__":
    main()