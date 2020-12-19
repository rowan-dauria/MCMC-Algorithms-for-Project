'''
This code helps illustrate how Hamiltonian Monte Carlo works.
This is not the full simulation.

The blue lines trace the red contours (contours along which the
the Hamiltonian(H) is conserved). This shows how the proposals
conserve H, which allows for distant proposals with a high
probability of acceptance. (RM Neal 2011,
www.arxiv.org/abs/1206.1901)


'''

import numpy as np
import matplotlib.pyplot as plt

#mass of simulated particle
m = 1.


# p is represented by a 1x3 array
#K = p^2/(2m)
def K(p):
    K = np.vdot(p,p)/ (2*m)
    return float(K)

#harmonic potential
def U(q):
    U = 1/2 * np.vdot(q,q)
    return float(U)

#Hamiltonian of the system
def H(q,p):
    return U(q) + K(p)

#dU/dq = q
def dUdq(q):
    return q


'''generates random 3d momentum vectoreach component is selected from
normal dist, where the stdev is sqrt(m)'''
def randp():
    randp = np. vstack(
        [np.random.normal(0, np.sqrt(m), 1),
        np.random.normal(0, np.sqrt(m), 1),
        np.random.normal(0, np.sqrt(m), 1)])
    return randp




'''qs and ps are appended to with each value of p and q from the
hmc updating trajectory '''
qs = np.vstack([0.,0.,0.])
ps = np.vstack([0.,0.,0.])

'''function to perform HMC algorithm, works for any number of dimensions,
as long as you change randp() and input vectors accordingly'''
'''this version has been altered to visualise the process,
usually the coords of the curved trajectories are not stored,
except the very last one.'''
def nDHMC(current_q, dt, nstep, qs, ps):
    q = current_q
    p = randp()
    current_p = p
    qs = np.append(qs,q, axis=1)
    ps = np.append(ps,p, axis=1)
    
    p = p - (dt/2) * dUdq(q)
    for i in range(1, nstep-1):
        q = q + dt * p/m
        qs = np.append(qs,q, axis=1)
        if i != nstep:
            p = p - dt * dUdq(q)
            ps = np.append(ps,p, axis=1)
    p = p - (dt/2) * dUdq(q)
    
    current_U = U(current_q)
    current_K = K(current_p)
    proposed_U = U(q)
    proposed_K = K(p)

#trying to conserve Hamiltonian. If it is conserved, there will be no change to
# U or K, so the exponent will be zero. 
    if np.random.uniform() < np.exp(current_U-proposed_U+current_K-proposed_K):
        return q, qs, ps
    else:
        print('FAIL')
        #delete trjactory of rejected proposal
        return current_q, qs[:,:(-nstep-1)], ps[:,:(-nstep-1)]


q = np.vstack([0,0,0])#starting point vector 
for i in range(3): #gives 3 hmc updates
    q, qs, ps = nDHMC(q, .071, 18, qs, ps)#change dt and nstep to see how
#they affect final output
    

qs, ps = qs[:,1:], ps[:,1:] #removing first pos/n and mom/m that was the used
#to create arrays and is not a simulated point, both np.vstack[0,0,0]
    

#Making data for H contour plots
x, y = np.linspace(-3.,3.,200), np.linspace(-3.,3.,200)
Hs = np.zeros((len(x), len(y)))
for i, xi in enumerate(x):
    for j, yj in enumerate(y):
        Hs[i][j] = H(xi,yj)


fig, axs = plt.subplots(1,3, figsize = (18,6), sharey=True)
axs[0].plot(qs[0,:], ps[0,:])
axs[0].contour(x,y,Hs,[.125,.25,.5,1,2,3,4],colors='red', alpha=.4)
axs[0].set_ylabel('p')
axs[0].set_xlabel('q')

axs[1].plot(qs[1,:], ps[1,:])
axs[1].contour(x,y,Hs,[.125,.25,.5,1,2,3,4],colors='red', alpha=.4)
axs[1].set_xlabel('q')

axs[2].plot(qs[2,:], ps[2,:])
axs[2].contour(x,y,Hs,[.125,.25,.5,1,2,3,4],colors='red', alpha=.4)
axs[2].set_xlabel('q')
plt.show()


