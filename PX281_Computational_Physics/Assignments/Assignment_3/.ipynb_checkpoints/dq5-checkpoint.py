# YOUR CODE HERE
''' Question 5'''
from math import sqrt, pi
import numpy as np
from scipy.stats import maxwell
import scipy.constants as pc
import matplotlib.pyplot as plt


def k(v):
    ''' Kinetic Energy'''
    return 1/2 * 4 * pc.u * (v**2)


def k_(En):
    ''' Reverse Kinetic Energy '''
    return np.sqrt(2*En / (4*pc.u))


def samples(T1, T2, mass):
    ''' Sampling Function '''
    a1 = np.sqrt((pc.k * T1)/mass)
    a2 = np.sqrt((pc.k * T2) / mass)
    return (maxwell.rvs(scale=a1, size=1000),
            maxwell.rvs(scale=a2, size=1000))


def doCollision(ncoll, sample1, sample2):
    ''' Collision Simulation '''
    for _ in range(ncoll):
        s1 = np.random.randint(len(sample1))
        s2 = np.random.randint(len(sample2))
        e_1 = k(sample1[s1])
        e_2 = k(sample2[s2])
        news = np.abs(e_1-e_2)/2
        if sample1[s1] > sample2[s2]:
            sample1[s1] = k_(e_1-news)
            sample2[s2] = k_(e_2+news)
        else:
            sample1[s1] = k_(e_1+news)
            sample2[s2] = k_(e_2-news)
    return np.append(sample1, sample2)


m = 4.0*pc.u
wh, ch = samples(290, 4.2, m)
final = doCollision(10000, wh, ch)
T = (np.mean(final)/2.0*sqrt(pi/2))**2 * m / pc.k
print(f"The mean temperature is {T}K")
plt.figure(figsize=(16, 9))
plt.title("Maxwell-Boltzmann Distribution after collision", fontsize=18)
plt.ylabel("Density", fontsize=15)
plt.xlabel("Speed in m/s", fontsize=15)
plt.hist(final, bins=50, density=True)
plt.grid(True)
plt.show()
