volTensorField gradU = fvc::grad(U) ;

//volSymmTensorField S = symm(gradU) ;

volTensorField S= gradU/2 + gradU.T()/2 ;

volTensorField R= gradU/2 - gradU.T()/2 ;

dimensionedTensor Itensor
( 
	"Identity", 
	dimensionSet(0, 0, 0, 0, 0, 0, 0), 
	tensor::I 
);

fvTensorMatrix AEqn
            (
               fvm::ddt(A)  + fvm::div(phi,A) - fvm::laplacian(nu_stab, A )
             ==
               ( A & S ) + ( S & A ) + ( A & R ) - ( R & A ) -  (1.0/lambda)* (  A - Itensor  )  - fvc::laplacian(nu_stab, A )
            );

AEqn.relax();
AEqn.solve();

// Eigen decomposition using a QR algorithm of Eigen library 
Eigen::Matrix3d D_eig(Eigen::Matrix3d::Zero(3,3));
Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;

forAll(A, cellI)
 {
   // Transfer D from OF to Eigen        
   const tensor& D_R = A[cellI];  
       
   D_eig(0,0)=D_R.xx();
   D_eig(1,1)=D_R.yy();
   D_eig(2,2)=D_R.zz();

   D_eig(0,1)=D_R.xy();
   D_eig(1,0)=D_R.xy();

   D_eig(0,2)=D_R.xz();
   D_eig(2,0)=D_R.xz();

   D_eig(1,2)=D_R.yz();
   D_eig(2,1)=D_R.yz();
    
   // Compute eigenvalues/vectors in Eigen

   eigSol.compute(D_eig);
   Eigen::Vector3d eival = eigSol.eigenvalues();
   Eigen::Matrix3d eivect = eigSol.eigenvectors();

   // Transfer eigenvalues/vectors from Eigen to OF 
   tensor& vecsR = eigVecs_[cellI];  

   vecsR.xx()=eivect(0,0);
   vecsR.yx()=eivect(1,0);
   vecsR.zx()=eivect(2,0);

   vecsR.xy()=eivect(0,1);
   vecsR.yy()=eivect(1,1);      
   vecsR.zy()=eivect(2,1);      

   vecsR.xz()=eivect(0,2);  
   vecsR.yz()=eivect(1,2);
   vecsR.zz()=eivect(2,2);
   
   if (eival(0) >= 0) {
   eigVals_[cellI].xx()=eival(0);
   } 
   else {
   eigVals_[cellI].xx()=0.;
   }

   if (eival(1) >= 0) {
   eigVals_[cellI].yy()=eival(1);
   } 
   else {
   eigVals_[cellI].yy()=0.;
   }

   if (eival(2) >= 0) {
   eigVals_[cellI].zz()=eival(2);
   } 
   else {
   eigVals_[cellI].zz()=0.;
   }

 }



A = eigVecs_ & eigVals_ & eigVecs_.T() ;


tau_p = (nu_p/lambda) * symm( A - Itensor);

tau_p.correctBoundaryConditions();




