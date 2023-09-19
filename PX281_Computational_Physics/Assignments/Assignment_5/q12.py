""" My code """
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import bisect
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt


""" required physical constants given in SI units"""
m_el = 9.1094e-31  # mass of electron in [kg]
hbar = 1.0546e-34  # Planck's constant over 2 pi [Js]
e_el = 1.6022e-19  # electron charge in [C]
L_bohr = 5.2918e-11  # Bohr radius [m]

# YOUR CODE HERE


def eqn(x, y, energy):
    """Schrodinger's Equation"""
    _ = x
    w_, dw_ = y
    f_0 = dw_
    f_1 = (-2.0 * m_el / (hbar * hbar)) * w_ * energy  # dw = -2m/h^2 Ew
    return [f_0, f_1]


def solve(energy, func):
    """Solving Function"""
    x_range = [0, L_bohr]
    x_eval = np.linspace(x_range[0], x_range[1], 1000)
    init = [0, 1]
    answer = solve_ivp(func, x_range, init, args=(energy,), t_eval=x_eval)
    return answer.y[0, -1]


gr_en = bisect(solve, 130 * e_el, 140 * e_el, args=(eqn,), xtol=1e-22) / e_el
print(f"The ground state energy is at {gr_en}eV")
