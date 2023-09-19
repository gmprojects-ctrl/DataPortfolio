import numpy as np
import scipy.integrate as sci 
import matplotlib.pyplot as plt
''''
Assuming K,my and q are  constant
The model is 

mu omega^2 w = - K d4/dx4 w  + q

Let M = mu omega^2

Now rewriting as 

w_1 = w
w_2 = w_1' = w''
w_3 = w_2' = w''
w_4 = w_3' = w''\'

Then we get coupled dif

w_1' = w_2
w_2' = w_3
w_3' = w_4
w_4' = - M/K w  + q/K

'''

def func(x,y,mu,K,q,omg):
    M = mu*omg*omg
    w1,w2,w3,w4 = y # Tuple unpacking
    return [w2,w3,w4,-(M/K)*w1 +q/K]
def main() -> None:
    S_0 = (0,0,0,0) # Inital Conditions 
    L=10
    N=100
    omg = 0.25
    mu = 1
    K = 2
    q = 1
    x = np.arange(0,L,1/N)
    t = np.arange(0,L,1/N)
    sol = sci.solve_ivp(func,(0,L),S_0,t_eval=x,args=[mu,K,q,omg])
    w   = sol.y[1]*np.sin(omg*t)
    fig = plt.figure()
    ax  = fig.add_subplot(121,projection='3d')
    ax.scatter(t,x,w,color='r')
    ax.set_xlabel("t",size=12)
    ax.set_ylabel("x",size=12)
    ax.set_zlabel("w(t,x)",size=12)
    ax2 = fig.add_subplot(122)
    ax2.scatter(x,sol.y[1])
    plt.show()


if __name__ =="__main__":
    main()