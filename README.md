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
Applying a kernel transformation as  and being  <img src="https://latex.codecogs.com/svg.image?\Theta&space;=&space;\ln(A)=U\ln(\Lambda)U^{T}" title="\Theta = \ln(A)=U\ln(\Lambda)U^{T}" />, the following equation defines the evolution of \Theta:

<img src="https://latex.codecogs.com/svg.image?\frac{D&space;\Theta&space;}{Dt}&space;=&space;\Omega&space;\Theta&space;-&space;\Theta&space;\Omega&space;&plus;2B&plus;&space;\frac{e^{-\Theta}&space;-&space;&space;I&space;}{\lambda}" title="\frac{D \Theta }{Dt} = \Omega \Theta - \Theta \Omega +2B+ \frac{e^{-\Theta} - I }{\lambda}" />

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
Here an eigen decomposition method is employed and the following approaches were followed:

   i) Numerically solving the transport equation for the eigen values
   
   a)Oldroyd-B
     
   b)FENE-P
     
   ii) Application of a pseudo-analitical solution
   
   a)Solution following the assumption: <img src="https://latex.codecogs.com/svg.image?\frac{D\Lambda&space;}{Dt}\approx&space;\frac{\partial&space;\Lambda}{\partial&space;t}" title="\frac{D\Lambda }{Dt}\approx \frac{\partial \Lambda}{\partial t}" />
solver:

   b) Above solution but with explicit addition of the <img src="https://latex.codecogs.com/svg.image?(u\cdot&space;\triangledown&space;)&space;(\Lambda&space;)" title="(u\cdot \triangledown ) (\Lambda )" />
     
## 5 - Direct prescription of the conformation matrix - A.
 
   i) Linear identity : 
    <img src="https://latex.codecogs.com/svg.image?A=2&space;\lambda&space;D" title="A=2 \lambda D" />
    For this identity 
    
   ii) Quadratic identity : 
    <img src="https://latex.codecogs.com/svg.image?A=2&space;\lambda&space;D&space;&plus;&space;(2\lambda)^{2}DD" title="A=2 \lambda D + (2\lambda)^{2}DD" />
    
   iii) Cubic identity : <img src="https://latex.codecogs.com/svg.image?A=2&space;\lambda&space;D&space;&plus;&space;(2\lambda)^{3}DDD" title="A=2 \lambda D + (2\lambda)^{3}DDD" />
   
   iv) Identity with convection: 
    <img src="https://latex.codecogs.com/svg.image?&space;A=2&space;\lambda&space;D&space;&plus;&space;(5\lambda)^{2}(u\cdot&space;\triangledown&space;)(D)" title=" A=2 \lambda D + (5\lambda)^{2}(u\cdot \triangledown )(D)" />

Tutorials are presented to test every solver. 


## Conclusions:
