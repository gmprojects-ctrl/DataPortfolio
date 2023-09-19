'''
Free fall simulation with measurement errors.

This example contains programming faults and also physics blunders!
Hint: Suppress gravity accelerations above 29 m/s^2.
'''

from random import gauss
import numpy as np
import matplotlib.pyplot as plt


def gravity(height, time):
    ''' Gravity Func'''
    grav = 2*height/time**2
    return grav if grav <= 29 else-1


def fallsim(attempts, height, heightError, time, timeError):
    ''' Fall Sim '''
    collector = np.array([])
    for _ in range(attempts):  # no counter variable needed
        distance = gauss(height, heightError)
        watch = gauss(time, timeError)
        collector = np.append(collector, gravity(distance, watch))
    return collector[collector != -1]


pisa = 58  # [m]
falltime = 3.4  # [s] in standard Earth gravity
herror = 0.5  # [m] uncertainty from where to drop to the floor
werror = 1.0  # [s] watches were rather uncertain
measurements = fallsim(10000, pisa, herror, falltime, werror)
plt.hist(measurements, 21)
plt.title('Gravity constant measurements')
plt.ylabel('measurements', fontsize=14)
plt.xlabel('g [ms$^{-2}$]', fontsize=14)
plt.show()
