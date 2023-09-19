# Solution to part (a)
# YOUR CODE HERE
""" My Code """
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# init = [x,dy,y,dy]


def func(t, y):
    """Func"""
    _ = t
    r = y[0] * y[0] + y[2] * y[2]
    f_0 = y[1]
    f_1 = -2.0 * (y[2] * y[2]) * y[0] * (1 - y[0] * y[0]) * np.exp(-r)
    f_2 = y[3]
    f_3 = -2.0 * (y[0] * y[0]) * y[2] * (1 - y[2] * y[2]) * np.exp(-r)
    return np.array([f_0, f_1, f_2, f_3])


def trajectory(impactpar, speed):
    """Trajectory"""
    init = [impactpar, 0.0, -2, speed]
    maxtime = 10 / speed
    time = np.linspace(0, maxtime, 300)
    ans = solve_ivp(func, [0, maxtime], init, t_eval=time)
    return (ans.y[0], ans.y[2])


x, y_ = trajectory(0.5, 0.25)
plt.figure(figsize=(16, 9))
plt.grid(True)
plt.title("y(t) against x(t)", fontsize=18)
plt.xlabel("Displacement in x", fontsize=18)
plt.ylabel("Displacement in y", fontsize=18)
plt.plot(x, y_)
plt.show()
