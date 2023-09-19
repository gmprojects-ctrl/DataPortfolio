import math
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt

# *********************************************
# original function: 
# \phi(x) = 2 \sin(3x) + 0.7 \sin(11x) 
# TODO 
# after addressing the other task below, 
# implementing another function for phi 
# *********************************************
def phifct(x):
    phi = 2.0*np.sin(3.0*x) + 0.7*np.sin(11.0*x)
    return phi

# setting up the problem 
N = 1025 # spatial discretisation, 
         # N is the number of mesh points for plotting 
xx = np.linspace(-math.pi,math.pi,N) 

# defining and plotting phi
phi = np.zeros(N)
for i in range(0,N):    
    phi[i] = phifct(xx[i]) 
fig = plt.figure()
plt.plot(xx,phi,'b')
plt.show()

# discrete Fourier transform
zphi = fft(phi) 
plt.plot(np.abs(zphi),'k')
plt.show()
# zphi is an array of length N, 
# \hat{\phi}(0) is the first entry in zphi(0), 
# then zphi(k) = \hat{\phi}(k) and 
# and  zphi(N-k) = \hat{\phi}(-k)

# cut off high frequencies, 
# in MATLAB, given a frequency n, 
# entries between n+1 and N-n+1 are set to zero 
# *********************************************
# TODO 
# increase n until the original function 
# is correctly recovered, ie the red curve is  
# on top of the blue curve in the last graph
# *********************************************
n = 11; 
zphiapprox = zphi;
for i in range(0,N):
    if (i>n and i<N-n):
        zphiapprox[i] = 0.0
plt.plot(np.abs(zphiapprox),'g')
plt.show()

# transform back and plot both phi and its approximation
phiapprox = ifft(zphiapprox) 
plt.plot(xx,phi,'b')
plt.plot(xx,phiapprox,'r')
plt.show()
