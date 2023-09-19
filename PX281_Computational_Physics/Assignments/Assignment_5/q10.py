# YOUR CODE HERE
""" My Code """
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def func(t, y, mu_):
    """My Func"""
    _ = t
    x, dx, y_, dy_ = y
    r = np.sqrt((x + mu_) ** 2 + y_ ** 2)
    s = np.sqrt((x - 1 + mu_) ** 2 + y_ ** 2)
    f_0 = dx
    f_1 = x + 2 * dy_ - (1 - mu_) * (x + mu_) / r ** 3
    f_1 = f_1 - mu_ * (x - 1 + mu_) / s ** 3
    f_2 = dy_
    f_3 = y_ - 2 * dx - (1 - mu_) * (y_) / r ** 3 - mu_ * y_ / s ** 3
    return np.array([f_0, f_1, f_2, f_3])


def satellite(init_, time_, mu_):
    """Satellite"""
    evalt = np.linspace(time_[0], time_[1], 200)
    return solve_ivp(func, time_, init_, args=(mu_,), t_eval=evalt)


mu = 0.01227471
init = np.array([0.994, 0.0, 0.0, -2.0015851])
time = [0, 18]
time_range = np.linspace(time[0], time[1], 200)
answer = satellite(init, time, mu)
# Plotting
plt.figure(figsize=(18, 12))
plt.subplot(131)
plt.grid(True)
plt.title("y(t) against x(t)", fontsize=18)
plt.xlabel("Displacement in x", fontsize=12)
plt.ylabel("Displacement in y", fontsize=12)
plt.plot(answer.y[0], answer.y[2], "r")
plt.subplot(132)
plt.grid(True)
plt.title("y(t) and x(t) against time", fontsize=18)
plt.xlabel("Time", fontsize=12)
plt.ylabel("Displacement", fontsize=12)
plt.plot(time_range, answer.y[1], "r", label="x(t)")
plt.plot(time_range, answer.y[2], "b", label="y(t)")
plt.legend(fontsize=12)
plt.subplot(133)
plt.grid(True)
plt.title("dy/dt against dx/t", fontsize=18)
plt.xlabel("Velocity in x", fontsize=12)
plt.ylabel("Velocity in y", fontsize=12)
plt.plot(answer.y[1], answer.y[3], "r")
plt.show()
