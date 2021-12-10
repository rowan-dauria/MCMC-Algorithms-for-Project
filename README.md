# MCMC-Algorithms-for-Project
These scripts form the basis of my final year university computing project. The simulation scripts use different type of Markov Chain Monte Carlo algorithms to simulate the wavefuction of a quantum particle in 3D. Essentially, both algorithms propose points where the particle might be found, and the statistical distributionn of these points is the simulated probablity distribution of the particle. 

The scripts are generalised to allow simulation of higher or lower dimension systems.

### nD_metropolis.py
Proposes possible positions for each particle using the Metropolis algorithm. This is a fundemental MCMC algorithm but it cannot make distant proposals, which makes it less effective for simulating mutli-peaked probablity distributions. 

### nD_HMC.py
Proposes points using the Hamiltonian Monte Carlo methdod. This is a more advanced MCMC algorithm that can propose more distant points, which is effective for simulating more complex probability distributions. However, there are more calculations required for each proposal. 


The Markov Chain Monte Carlo algorithms used in my final year computing project. 

You will need matplotlib and numpy to run these scripts. 
