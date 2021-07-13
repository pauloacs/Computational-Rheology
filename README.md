# Computational-Rheology

Here some solvers in OpenFOAM able to simulate viscoelastic fluids are presented. The solvers are divided in 2 batches:

1 - Evolution of the conformation matrix 
Here the Oldroyd-B constitutive model was used. 

a) Direct application of the differential equation for the conformation matrix:

<img src="https://latex.codecogs.com/svg.image?&space;\frac{DA}{Dt}&space;=&space;AS&space;&plus;&space;SA&space;&plus;&space;AR&space;-RA&space;-&space;\frac{A&space;-&space;I&space;}{\lambda}&space;" title=" \frac{DA}{Dt} = AS + SA + AR -RA - \frac{A - I }{\lambda} " />

solver: 
b) With kernel transformation

<img src="https://latex.codecogs.com/svg.image?&space;\frac{D&space;\Theta&space;}{Dt}&space;=&space;\Omega&space;&space;\Theta&space;-&space;\Theta&space;&space;\Omega&space;&plus;2B&plus;&space;\frac{&space;e^{-\Theta}&space;-&space;I&space;&space;}{\lambda}&space;" title=" \frac{D \Theta }{Dt} = \Omega \Theta - \Theta \Omega +2B+ \frac{ e^{-\Theta} - I }{\lambda} " />

solver:

2- Evolution of the eigen values and eigen vectors of the conformation matrix - Here various considerations have been tested.

1- Collins'_method_exp - Method from "Numerical approach to simulating turbulent flow of a viscoelastic polymer solution" by T. Vaithianathan and Lance R. Collins with exponential kernel transformation

2- Collins' approach

3- Assuming A=lambda*D with further mathematical simplifications

4- Using eigen decompotion on D (deformation rate tensor) to use the eigen vectors as the principal directions

  i - Numerically solving principal deformations dif. equation
  
  ii - Approx analytical solution to the evolution of principal deformations
  
    a)  ..
    b) incremental

6 - 

Tutorials are presented to test every solver. 
