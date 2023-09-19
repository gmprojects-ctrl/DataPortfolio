import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("output.txt")
time = np.unique(data[:,0])
print()
plt.xlabel("x",size=24)
plt.ylabel("Y(x,t)",size=24)
for i in time:
    x = data[data[:,0]==i][:,1]
    y = data[data[:,0]==i][:,2]
    plt.plot(x,y,label=f"Time = {i}")
plt.legend(fontsize=12)
plt.show()
