# Computational-Rheology

Here some solvers in OpenFOAM able to simulate viscoelastic fluids are presented. This implementations were done in v7 of OpenFOAM. In order to use some convective schemes and boundary conditions not typically present in openFOAM, in each solver in the file Make/options some files from ... are included. 
It is possible to compile most of the solvers(except the ones after II-3) by removing this, but the convective schemes in system/fvSchemes and the boundary conditions for tau_p must be changed in the tutoriais, but in order for the simulations to give the best results, rheoTool must be installed in openfoam 7. 

The solvers are divided in 2 batches:

# I - Evolution of the conformation matrix 
Here the *Oldroyd-B constitutive model* was used. 

## a) Direct application of the differential equation for the conformation matrix:
The spatial and temporal evolution of the conformation matrix can be given by the following differential equation:
<img src="https://latex.codecogs.com/svg.image?&space;\frac{DA}{Dt}&space;=&space;AS&space;&plus;&space;SA&space;&plus;&space;AR&space;-RA&space;-&space;\frac{A&space;-&space;I&space;}{\lambda}&space;" title=" \frac{DA}{Dt} = AS + SA + AR -RA - \frac{A - I }{\lambda} " />

solver: 
## b) With kernel transformation
Applying a kernel transformation as  and being <img src="https://latex.codecogs.com/svg.image?\Theta&space;=&space;\ln(A)=U\ln(\Lambda)U^{T}" title="\Theta = \ln(A)=U\ln(\Lambda)U^{T}" />, the following equation defines the evolution of \Theta.
<img src="https://latex.codecogs.com/svg.image?&space;\frac{D&space;\Theta&space;}{Dt}&space;=&space;\Omega&space;&space;\Theta&space;-&space;\Theta&space;&space;\Omega&space;&plus;2B&plus;&space;\frac{&space;e^{-\Theta}&space;-&space;I&space;&space;}{\lambda}&space;" title=" \frac{D \Theta }{Dt} = \Omega \Theta - \Theta \Omega +2B+ \frac{ e^{-\Theta} - I }{\lambda} " />

solver:

# II- Evolution of the eigen values and eigen vectors of the conformation matrix - 
This whole approach was encouraged by T. Vaithianathan and Lance R. Collins in "Numerical approach to simulating turbulent flow of a viscoelastic polymer solution". The proposed implementations such as various other implementations were considered. 

The conformation matriz - A - is a 

## 1 - Pure implementation
solver:
## 2 - Kernel transformation application to the above algorithm
solver:
## 3 - Evolution of each component of the eigen vectors
solver:
## 4 - Eigen Vectores extracted from deformation rate tensor
solver:
## 5 - Direct presciption of the conformation matrix - A.

Tutorials are presented to test every solver. 
