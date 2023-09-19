# Solution to part (b)
# YOUR CODE HER
''' My Code'''
import math


def scatterangles(allb, speed):
    ''' Scatter Angles '''
    angles = np.array([])
    for i in allb:
        init = [i, 0.0, -2, speed]
        maxtime = 10 / speed
        time = np.linspace(0, maxtime, 300)
        ans = solve_ivp(func, [0, maxtime], init, t_eval=time)
        angle = math.atan2(ans.y[3][-1], ans.y[1][-1])
        angles = np.append(angles, angle)
    return angles


impact_range = np.arange(-0.2, 0.2, 0.001)
angles_ = scatterangles(impact_range, 0.1)
plt.figure(figsize=(16, 9))
plt.grid(True)
plt.title("Scatterangles against Impact Range", fontsize=18)
plt.xlabel("Impact Range", fontsize=18)
plt.ylabel("Scatter Angle", fontsize=18)
plt.scatter(impact_range, angles_)
plt.show()
