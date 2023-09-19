import math 
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import default_rng
from scipy import sparse 
from scipy.sparse.linalg import spsolve 

def generate_mesh (lx1, lx2, n):
    # problem dependent generation of a mesh in 2D
    # as part of the rectangle R = [lx1(1) lx1(2)] \times [lx2(1) lx2(2)]
    # with grid parameters hx, hy 
    # given by n = (n[0],n[1]) points in the coordinate directions 
    # 
    # input data:
    #   lx1 : interval in x1 direction 
    #   lx2 : interval in x2 direction 
    #   n   : n[0] n[1], number of mesh points in x1 and x2 direction 
    # 
    # output data:
    #   X1(0,...,n[0]-1) : vector of the x1 coordinates of the vertices
    #   X2(0,...,n[1]-1) : vector of the x2 coordinates of the vertices
    #   M  : mesh, n[0] x n[1] array with vertex number for each vertex 
    #   ct : number of interior vertices 
    # 
    # B Stinner 2020

    X1 = np.linspace(lx1[0],lx1[1],n[0])
    X2 = np.linspace(lx2[0],lx2[1],n[1])
    M = np.zeros((n[0],n[1]), dtype=int) 
    ct = 0
    for i in range(0,n[0]): 
        for j in range(0,n[1]): 
            flg = 1 # taking the whole domain 
            if (flg > 0):  # (X1(i),X2(j)) is a mesh point 
                M[i,j] = ct
                ct = ct+1
            else: 
                M[i,j] = flg

    return X1, X2, M, ct

def generate_mesh_L (lx1, lx2, n): 
    # problem dependent generation of a mesh in 2D
    # as part of the rectangle R = [lx1(1) lx1(2)] \times [lx2(1) lx2(2)]
    # with grid parameters hx, hy 
    # given by n = [n(1) n(2)] points in the coordinate directions 
    # 
    # input data:
    #   lx1 : interval in x1 direction 
    #   lx2 : interval in x2 direction 
    #   n   : n[0] n[1], number of mesh points in x1 and x2 direction 
    # 
    # output data:
    #   X1(1,...,n[0]) : vector of the x1 coordinates of the vertices
    #   X2(1,...,n[1]) : vector of the x2 coordinates of the vertices
    #   M  : mesh, n[0] x n[1] array with vertex number for each vertex 
    #   ct : number of interior vertices 
    # 
    # B Stinner 2020
    
    X1 = np.linspace(lx1[0],lx1[1],n[0])
    X2 = np.linspace(lx2[0],lx2[1],n[1])
    M = np.zeros((n[0],n[1]), dtype=int) 
    ct = 0
    tol = 1.0e-5
    for i in range(0,n[0]): 
        for j in range(0,n[1]): 
            if (X1[i]>0.6*math.pi+tol and X2[j]>0.35*math.pi+tol): 
                flg = -1
            else:
                flg = 1 
            if (flg > 0):  # (X1(i),X2(j)) is a mesh point 
                M[i,j] = ct
                ct = ct+1
            else: 
                M[i,j] = flg
    
    return X1, X2, M, ct

def init_u (X1,X2,M,di): 
    # problem dependent function for initial  values
    # 
    # loop over the mesh points to unitialise u 
    # 
    # input data:
    #   X1, X2: coordinates of the mesh points
    #   M: mesh, n[0]xn[1] array with vertex number for each vertex in the
    #      domain, contains a nonpositive number if the point it not in the domain
    #  di: number of vertices 
    # 
    # output data:
    #   Uval: initial values for u in the vertices.
    # 
    ub = np.array([42.8576,0.0184])  
    uval = np.zeros((di,2)) 
    rng = default_rng() 
    n = M.shape  
    ct = 0 
    for i in range(0,n[0]): 
        for j in range(0,n[1]): 
            if (M[i,j] >= 0):  
                uval[ct,0] = ub[0] * (1 + 0.05*np.cos(50*X1[i])*np.cos(50*X2[j])) * (0.97 + 0.06*rng.uniform())  
                uval[ct,1] = ub[1] * (0.97 + 0.06*rng.uniform()) 
                ct = ct+1 
    return uval

def react_fct (uval):
    # problem dependent function for reaction terms
    # 
    # **************** TODO ******************
    # implement the correct reaction function
    # ****************************************
    # 
    alpha = 1500 # Values from 3b
    rb = 0.015
    mu = 35
    u1 = uval[0]
    u2 = uval[1]
    x = rb*(1+(u1*u1)/u2) - mu*u1
    y = rb*u1*u1 - alpha*u2
    return x,y

def vis_fct (fig, XX1, XX2, U1, U2): 
    # visualisation of the solution 
    plt.clf() 
    axs1 = fig.add_subplot(2, 2, 1, projection='3d')
    axs2 = fig.add_subplot(2, 2, 2, projection='3d')
    axs3 = fig.add_subplot(2, 2, 3)
    axs4 = fig.add_subplot(2, 2, 4)
    axs1.plot_surface(XX1, XX2, U1, cmap='viridis', edgecolor='none')
    axs1.set_title('u_1 (side view)')
    axs2.plot_surface(XX1, XX2, U2, cmap='viridis', edgecolor='none')
    axs2.set_title('u_2 (side view)')
    axs3.pcolor(XX1, XX2, U1, cmap='viridis')
    axs3.set_title('u_1 (top view)')
    axs4.pcolor(XX1, XX2, U2, cmap='viridis')
    axs4.set_title('u_2 (top view)')
    plt.draw()
    plt.pause(0.001)
    return

def rds_solve (tau, no_steps, lx1, lx2, n, dtvis):
    # solver for a reaction diffusion system of the form 
    #   d_t u(x,t) - D(x,t) \Laplacian u(x,t) + r(u(x,t)) = s(x,t) 
    # in 2D with Neumann boundary condition 
    # using finite difference techniques 
    # 
    # input data:
    #   tau, no_steps          : time step, number of time steps
    #   lx1 = [lx1(1),lx1(2)]  : interval in x1 direction 
    #   lx2 = [lx2(1),lx2(2)]  : interval in x2 direction 
    #   n = [n[0],n[1]]        : number of grid points in x1 and x2 direction
    #   dtvis                  : visualisation time step (switched off if <0) 
    # 
    # output data: 
    #   X1, X2, M: mesh
    #   uold: solution values 
    # 
    # example call: 
    # rds_solve (2.0e-3, 300, [0 pi], [0 pi], [80 80], 1.0e-2)
    # 
    # B Stinner 2020
    
    # problem parameters 
    d1 = 0.3 
    d2 = 90.0
    ub1 = 42.8576 
    ub2 = 0.0184 
    
    # visualisation flags 
    if dtvis>0:   # solution is visualised
        vis_flag = 1 
    else:
        vis_flag = 0 
    
    # initialise mesh and u values 
    # 
    # **************** TODO ******************** 
    # use generate_mesh_L for the last question 
    # ****************************************** 
    # 
    X1, X2, M, di = generate_mesh (lx1, lx2, n) 
    uold = init_u (X1, X2, M, di) 
    
    # initial visualisation
    if (vis_flag == 1): 
        fig = plt.figure()
        # fig, axs = plt.subplots(2, 2)
        U1 = np.zeros((n[0],n[1])) 
        U2 = np.zeros((n[0],n[1])) 
        for i in range(0,n[0]): 
            for j in range(0,n[1]):
                ent = M[i,j] 
                if (ent >= 0):
                    U1[i,j] = uold[ent,0] 
                    U2[i,j] = uold[ent,1] 
                else: # outside of the domain, use dummy value
                    U1[i,j] = ub1 
                    U2[i,j] = ub2 
        XX1, XX2 = np.meshgrid(X1,X2) 
        vis_fct (fig, XX1, XX2, U1, U2)
        timevis = dtvis
    
    # assembly of the stiffness matrix 
    A = sparse.lil_matrix((di,di)) 
    hx1 = (lx1[1]-lx1[0])/(n[0]-1) 
    hx2 = (lx2[1]-lx2[0])/(n[1]-1) 
    for i in range(0,n[0]): 
        for j in range(0,n[1]): 
            ent = M[i,j]
            if (ent >= 0): # we have a linear equation in this point 
                # diagonal entry
                A[ent,ent] = A[ent,ent] + 2.0 * (1.0/(hx1*hx1) + 1.0/(hx2*hx2))
                # off-diagonal, north 
                if (j == n[1]-1): # Neumann boundary vertex 
                    ent3 = M[i,j-1] 
                    A[ent,ent3] = A[ent,ent3] - 1.0/(hx2*hx2)
                else:
                    ent2 = M[i,j+1] 
                    if (ent2 < 0): # Neumann boundary vertex 
                        ent3 = M[i,j-1] 
                        A[ent,ent3] = A[ent,ent3] - 1.0/(hx2*hx2)
                    else:
                        A[ent,ent2] = A[ent,ent2] - 1.0/(hx2*hx2)
                # off-diagonal, south
                if (j == 0): # Neumann boundary vertex 
                    ent3 = M[i,j+1] 
                    A[ent,ent3] = A[ent,ent3] - 1.0/(hx2*hx2)
                else:
                    ent2 = M[i,j-1] 
                    if (ent2 < 0): # Neumann boundary vertex 
                        ent3 = M[i,j+1] 
                        A[ent,ent3] = A[ent,ent3] - 1.0/(hx2*hx2)
                    else:
                        A[ent,ent2] = A[ent,ent2] - 1.0/(hx2*hx2)
                # off-diagonal, west
                if (i == 0): # Neumann boundary vertex 
                    ent3 = M[i+1,j] 
                    A[ent,ent3] = A[ent,ent3] - 1.0/(hx1*hx1)
                else:
                    ent2 = M[i-1,j] 
                    if (ent2 < 0): # Neumann boundary vertex 
                        ent3 = M[i+1,j] 
                        A[ent,ent3] = A[ent,ent3] - 1.0/(hx1*hx1)
                    else:
                        A[ent,ent2] = A[ent,ent2] - 1.0/(hx1*hx1)
                # off-diagonal, east
                if (i == n[0]-1): # Neumann boundary vertex 
                    ent3 = M[i-1,j] 
                    A[ent,ent3] = A[ent,ent3] - 1.0/(hx1*hx1)
                else:
                    ent2 = M[i+1,j] 
                    if (ent2 < 0): # Neumann boundary vertex 
                        ent3 = M[i-1,j] 
                        A[ent,ent3] = A[ent,ent3] - 1.0/(hx1*hx1)
                    else:
                        A[ent,ent2] = A[ent,ent2] - 1.0/(hx1*hx1)
    A = A.tocsr()
    # other storage 
    D = sparse.eye(di)
    b = np.zeros((di,2))
    unew = np.zeros((di,2)) 
    
    # time loop 
    for m in range (1,no_steps+1):
        time = m*tau
        # assemble, reaction, semi-implicit
        for ct in range(0,di):
            rhs0, rhs1 = react_fct (uold[ct,:]) 
            b[ct,0] = uold[ct,0] + tau * rhs0 
            b[ct,1] = uold[ct,1] + tau * rhs1 
        # solve SLE for new u value 
        x = spsolve(A, b)
        unew[:,0] = spsolve((D + d1*tau*A),b[:,0]) 
        unew[:,1] = spsolve((D + d2*tau*A),b[:,1]) 
        # copy over    
        uold = unew 
        # visualise if desired
        if time >= timevis: 
            if (vis_flag == 1): 
                for i in range(0,n[0]): 
                    for j in range(0,n[1]):
                        ent = M[i,j] 
                        if (ent >= 0):
                            U1[i,j] = uold[ent,0] 
                            U2[i,j] = uold[ent,1] 
                        else: # outside of the domain, use dummy value
                            U1[i,j] = ub1 
                            U2[i,j] = ub2 
                vis_fct (fig, XX1, XX2, U1, U2)
            timevis = timevis + dtvis
        print('After step ', m, ', time: ', time)
    
    if (vis_flag == 1): # final visualisation 
        for i in range(0,n[0]): 
            for j in range(0,n[1]):
                ent = M[i,j] 
                if (ent >= 0):
                    U1[i,j] = uold[ent,0] 
                    U2[i,j] = uold[ent,1] 
                else: # outside of the domain, use dummy value
                    U1[i,j] = ub1 
                    U2[i,j] = ub2 
        vis_fct (fig, XX1, XX2, U1, U2)
        plt.show() 
    
    return X1, X2, M, uold 

# **************** TODO ******************** 
# run rds_solve with the correct parameters 
# ****************************************** 
# 
tau = 0.001
no_steps = 2500
lx1 = np.array([0,math.pi])
lx2 = np.array([0,math.pi])
n = np.array([80,80])
dtvis = 0.1
X1, X2, M, U1 = rds_solve (tau, no_steps, lx1, lx2, n, dtvis) 
