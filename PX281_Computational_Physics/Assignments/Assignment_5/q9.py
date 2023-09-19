# YOUR CODE HERE
""" My Code """
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def rateq(t, y, k_1, k_2):
    """Our literal ODEs"""
    _ = t
    dy_1 = -k_1 * y[0]
    dy_2 = k_1 * y[0] - k_2 * y[1]
    return [dy_1, dy_2]


def rateEqns(initial, time, k_1, k_2):
    """Rate Eqn"""
    return solve_ivp(rateq, time, initial, args=(k_1, k_2))


# Actual Data
k1 = 0.2
k2 = 0.8
y10 = 100
y20 = 0
init = np.array([y10, y20])
answer = rateEqns(init, [0, 20], k1, k2)
loss = y10 - answer.y[0] - answer.y[1]
plt.figure(figsize=(16, 9))
plt.grid(True)
plt.title("y1 & y2 and Loss feed against time", fontsize=18)
plt.xlabel("Time", fontsize=18)
plt.ylabel("Displacement", fontsize=18)
plt.plot(answer.t, answer.y[0], "r+", label="y1")
plt.plot(answer.t, answer.y[1], "bx", label="y2")
plt.legend(fontsize=18)
plt.figure(figsize=(16, 9))
plt.grid(True)
plt.title(" Loss feed against time", fontsize=18)
plt.xlabel("Time", fontsize=18)
plt.ylabel("Loss", fontsize=18)
plt.plot(answer.t, loss, "go", label="Loss")
plt.legend(fontsize=18)
plt.show()
