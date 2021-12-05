import numpy as np
import matplotlib.pyplot as plt

'''
Code to (roughly) simulate wave function of a quantum particle trapped in a
3D harmonic potential well with mass m, position q and momentum p using
the Hamiltonian algorithm.

Overall Method:
1. Define the field interacting with the simulated particle. 
2. Find its associated potential (in this case harmonic potential).
3. Choose a starting position for particle.
4. Propose a new position using the Hamiltonian algorithm.
5. Accept or reject the new position according to the energy associated with it. 
    Lower energy = accept (move to new position), higher energy = reject (particle doesn't move).
6. Record accepted points.
7. Plot a histogram of the points. 

Algorithm Method:
1. Particle is in position A (current_q).
2. Generate a random momentum for the particle (current_p) using a probability distribution defined by
    Hamiltonian algorithm. When position and momentum are defined, particle has energy current_E
3. Use "leapfrog" method to iterate to a proposed position (q) and momentum (p) that has energy 
    E such that E = current_E. In practise E ~ current_E, as leapfrog steps => 0, current_E = E. 
4. If E <= current_E, the position q is accepted and recorded. If E > current_E, there is a finite
        probabilty it will be accepted, but probably will be rejected and discarded. p is discarded. 
5. If q is accepted, repeat with current_q = q, otherwise current_q is unchanged. 
'''

m = 1 #value is 1 for simplicity

# p is represented by a vertical matrix, gives p^T * p as a scalar
def K(p): # K(p) is the kinetic energy of the particle
    K = 1/2 * np.vdot(p,p)/ m
    return float(K)


def U(q): #U(q) is the potential energy of the particle
    U = 1/2 * np.vdot(q,q)
    return float(U)

#qT * q is essentially x^2 + y^2 + ... so we can do normal partial differentials on it
def dUdq(q):
    return q


def randp(): #randp() proposes a new particle momentum, and a new particle position is chosen accordingly to keep the total energt (U+K) constant
    randp = np. vstack(
        [np.random.normal(0, np.sqrt(m), 1),
        np.random.normal(0, np.sqrt(m), 1),
        np.random.normal(0, np.sqrt(m), 1)])
    return randp

#the Hamiltonian algorithm      
def nDHMC(current_q, dt, nstep):
    q = current_q
    p = randp()
    current_p = p
    
    p = p - (dt/2) * dUdq(q)
    for i in range(1, nstep-1):
        q = q + dt * p/m
        if i != nstep:
            p = p - dt * dUdq(q)
    p = p - (dt/2) * dUdq(q)
    
    current_U = U(current_q)
    current_K = K(current_p)
    proposed_U = U(q)
    proposed_K = K(p)

    if np.random.uniform() < np.exp(current_U-proposed_U+current_K-proposed_K):
        return q
    else:

        return current_q

qs = np.vstack([0,0,0])
q=qs

for a in range(50000):
    q = nDHMC(q, .071, 18)
    qs=np.append(qs,q, axis = 1) # axis=1 -> 3d , axis=0 -> 1d

# analytic() defines the analytic (predicted by theory) function we are trying to model
def analytic(x):   # define analytic solution
    return (np.exp(-(x**2))) / (np.pi**0.5)



###for 3d plotting### -- the plotting dimensions were chosen arbitrarily
xdata = np.array(qs[0,:]).ravel()
ydata = np.array(qs[1,:]).ravel()
zdata = np.array(qs[2,:]).ravel()

fig, axs = plt.subplots(1,3, sharey = True, figsize = (9,3))

axs[0].hist(xdata, np.linspace(-3,3,100), density = True)
axs[0].plot(np.linspace(-3,3,100),analytic(np.linspace(-3,3,100)))

axs[1].hist(ydata, np.linspace(-3,3,100), density = True)
axs[1].plot(np.linspace(-3,3,100),analytic(np.linspace(-3,3,100)))

axs[2].hist(zdata, np.linspace(-3,3,100), density = True)
axs[2].plot(np.linspace(-3,3,100),analytic(np.linspace(-3,3,100)))

plt.show()
