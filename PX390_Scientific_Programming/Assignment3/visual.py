import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("output.txt")

time = data[:,0] 
x    = data[:,1]
y    = data[:,2]

fig = plt.figure()
ax  = fig.add_subplot(projection='3d')
ax.scatter(time,x,y,color='r')
ax.set_xlabel("t",size=12)
ax.set_ylabel("x",size=12)
ax.set_zlabel("Y(t,x)",size=12)
plt.show()
