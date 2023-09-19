'''Assignment 3'''
# function first, then main script for plotting
# YOUR CODE HERE
import numpy as np
import matplotlib.pyplot as plt
def fourier(xt,n):
    '''Fourier Function in Action'''
    count = 0
    s = lambda x,i: 1/(2*i-1)*np.sin((2*i-1)*2*np.pi*x)
    for i in range(1,n+1):
        count+=s(xt,i)
    return count*(4/np.pi)


points = np.linspace(0,1,200)
N = (3,11,40)
plt.figure(figsize=(16,9))
plt.title("Fourier Approximation of Step Function against x/T ", fontsize=24)
plt.xlabel("x/T",fontsize=18)
plt.ylabel("Approximated values of the Step Function",fontsize=18)
plt.grid(True)
for j in N:
    plt.plot(points,fourier(points,j),label=f"Fourier Approximation at n={j}")
plt.legend(fontsize=18)
plt.show()
