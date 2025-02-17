# ClassicalFluidSimulation

Classical Fluid Simulation for 1D and 2D analysis

The analyzed physical system is the classical fluid system. In here, the models of rigid rods and rigid disks will be used to represent the 1D and 2D systems of classical fluids. 
In order to simulate this system, I used Python. 
The tecnique used was the Monte Carlo Steps method to evolve the system and make it behave like a classical fluid. 
Also, very important for the simulation, is the establishment of periodic boundary conditions in all systems, both 1D as 2D, so that we simulate an “infinite” fluid. This is possible by making the boundary condition of a side being the state of the opposite side.

The parameters were pre-established by the system I wanted to simulate, following metrics that would make the calculation more coherent.
### Parameters:
* G(x) is the radial distribution function
* L, Lx and Ly are the size parameters (1D and 2D)
* N is the ammount of particles
* σ is the radius os the particles

ATTENTION: Due to difficulties previous to achieving the simulations, there are some commented parts in the Monte Carlo step for the 2D model that I did not refine after I had the desired results. But these parts are redundant and are not supposed to be uncommented.


If you have any problem or doubt during the local execution of some of those files, you can always contact me on one of my emails:

Gmail: [geraldo.siqueiraneto2000@gmail.com](mailto:geraldo.siqueiraneto2000@gmail.com)

Outlook: [geraldo.siqueiraneto@hotmail.com](mailto:geraldo.siqueiraneto@hotmail.com)
