# Computational-Rheology

Here some solvers in OpenFOAM able to simulate viscoelastic fluids are presented. The solvers are divided in 2 batches:

1 - Evolution of the conformation matrix

a) Direct application of the differential equation for the conformation matrix:


<img src="https://bit.ly/3i7NgFw" align="center" border="0" alt=" \frac{D \Theta }{Dt} = \Omega  \Theta - \Theta  \Omega +2B+ \frac{ e^{-\Theta}-I  }{\lambda} " width="249" height="46" />

solver: 
b) With kernel transformation
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

Tutorials:
