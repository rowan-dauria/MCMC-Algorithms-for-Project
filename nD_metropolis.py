import numpy as np
import matplotlib.pyplot as plt
'''It is advisable to use something other than np matrix'''
'''Vectors in this code are np.vstacks'''
x0 = -3
xn = 3
npoints = 60
nsteps = float(npoints - 1)
lattice = np.linspace(x0,xn,npoints)
spacing =(xn-(x0))/ nsteps
m=1
timestep = 100/60

def V(qa, qb):
    q = (qa+qb)/2
    V = 1/2 * np.vdot(q,q)
    return float(V)

def V1(qa, qb):
    q = (qa+qb)/2
    V = 1/2 * np.vdot(q,q) + q[0] - q[1]
    return float(V)
q = np.array([1.,1.,1.])
print(V1(q,q))

def Action(path, dt):
    Action = 0
    for i in range(1,path.shape[1]):
        q0, q1 = path[:,i-1], path[:,i]
        dq = q1-q0
        Ek =float(.5 * np.vdot(dq,dq) / (dt**2))
        Action += dt*(Ek + V1(q0, q1))
    return Action

def metropolis_sweep(path):
    p = .85 #perturbation size
    for i in range(path.shape[1]): #iterates along length of path
        new_path = path.copy()
        new_path[:,i] = path[:,i] + 2.*p * np.random.rand(1,path.shape[0]) - p
        if i == 0:
            dS = Action(new_path[:,i:i+2],timestep) - Action(path[:,i:i+2],timestep)
        elif i == path.shape[1]-1:
            dS = Action(new_path[:,i-1:i+1],timestep) - Action(path[:,i-1:i+1],timestep)
        else:
            dS = Action(new_path[:,i-1:i+2],timestep) - Action(path[:,i-1:i+2],timestep)
        #dS = Action(new_path, timestep) - Action(path, timestep)
        if dS < 0 or np.exp(-dS) > np.random.uniform(0,1): # is th second cond. correct?
            path[:,i] = new_path[:,i]#update the path
        else:
             path[:,i]=path[:,i]#keep orig path
    return path
        
###This code takes the same amount of time to run in 3d
###as it does in 1d!!
path = np.zeros((3,npoints))
print(metropolis_sweep(path))

for i in range(500):
    path = metropolis_sweep(path)

Ncf = 10000
Nrep = 10
paths = path.copy()

for i in range(1,Nrep+1):
    paths = path.copy()
    filename = 'x&y60pt-dt=1.67' + str([i]) + '.csv'

    for i in range(Ncf):
        path = metropolis_sweep(path)
        paths = np.concatenate((paths, path), axis=1)
    print(paths.T)
    print(filename)
    np.savetxt(filename, paths.T, delimiter = ',') 




def analytic(x):   # define analytic solution
    top = np.exp(-(x**2)/2)
    bottom = np.pi**0.25
    y = (top/bottom)**2
    return y





















    


