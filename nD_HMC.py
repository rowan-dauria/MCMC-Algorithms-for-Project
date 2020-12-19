import numpy as np
import matplotlib.pyplot as plt

'''Code to (roughly) simulate wave function of a particle trapped in a
3D harmonic potential well'''

m = 1 #value is 1 for simplicity



# p is represented by a vertical matrix, gives p^T * p as a scalar
def K(p):
    K = 1/2 * np.vdot(p,p)/ m
    return float(K)


def U(q):
    U = 1/2 * np.vdot(q,q)
    return float(U)

#qT * q is essentially x^2 + y^2 + ... so we can
#do normal partial diffs on it
def dUdq(q):
    return q

#q_test = np.matrix([2,3]).T
#p=p-.3*dUdq(q_test)
def randp():
    randp = np. vstack(
        [np.random.normal(0, np.sqrt(m), 1),
        np.random.normal(0, np.sqrt(m), 1),
        np.random.normal(0, np.sqrt(m), 1)])
    return randp


#t = randp() - .3 * dUdq(q_test)
#print(t)
     
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

q0 = np.vstack([0,0,0])
qs=q0
q=q0



for a in range(50000):
    q = nDHMC(q, .071, 18)
    qs=np.append(qs,q, axis = 1) # axis=1 -> 3d , axis=1 -> 1d

'''defines the analytic function we are trying to model'''
def analytic(x):   # define analytic solution
    return (np.exp(-(x**2))) / (np.pi**0.5)



###for 3d plotting### -- xyz were chosen arbitrarily
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
