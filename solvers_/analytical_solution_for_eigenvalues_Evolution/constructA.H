volTensorField gradU = fvc::grad(U) ;
volSymmTensorField D = symm(gradU) ;

scalar t = runTime.value(); 

// Eigen decomposition using a QR algorithm of Eigen library 
Eigen::Matrix3d D_eig(Eigen::Matrix3d::Zero(3,3));
Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;

forAll(D, cellI)
 {
   // Transfer D from OF to Eigen        
   const symmTensor& D_R = D[cellI];  
       
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
   
   if (t==runTime.startTime().value())
   {
	eigVals_[cellI].xx()=eival(0);
   	eigVals_[cellI].yy()=eival(1);
   	eigVals_[cellI].zz()=eival(2);
   }

 }


/*
 // Eigen decomposition using the iterative jacobi algorithm 

  forAll(D, cellI)
    {
       int N=3;
       int NROT=0;
       jacobi(D[cellI], N, eigVals_[cellI], eigVecs_[cellI], NROT);

    }
*/

/*
if (t==runTime.startTime().value())
{
	volScalarField f_r = (L2_ -3)/(L2_ - tr(eigVecs_ & eigVals_ & eigVecs_.T()) ) ;
}
*/

volScalarField f_r = (L2_ -3)/(L2_ - tr(eigVecs_ & eigVals_ & eigVecs_.T()) ) ;

volTensorField B= ( eigVecs_.T() & D & eigVecs_ );

//Evolute vals in time with analalytical solution for the transport of eigVals



const scalar eig1_0 = 0;
const scalar eig2_0 = 0;
const scalar eig3_0 = 0;


forAll(B, cellI)
 {
	tensor& B_ = B[cellI];
	scalar& fr = f_r[cellI];  	

	eigVals_[cellI] *= 0.;

	eigVals_[cellI].xx() = 1.0/(fr - 2.0 * B_.xx() * lambda.value()) + (eig1_0 - 1.0/(fr - 2.0 * B_.xx() * lambda.value() ))* Foam::exp( (2.0 * B_.xx() * lambda.value() - fr ) * t/lambda.value() );

	eigVals_[cellI].yy() = 1.0/(fr - 2.0 * B_.yy() * lambda.value()) + (eig2_0 - 1.0/(fr - 2.0 * B_.yy() * lambda.value() ))* Foam::exp( (2.0 * B_.yy() * lambda.value() - fr ) * t/lambda.value() );

	eigVals_[cellI].zz() = 1.0/(fr - 2.0 * B_.zz() * lambda.value()) + (eig3_0 - 1.0/(fr - 2.0 * B_.zz() * lambda.value() ))* Foam::exp( (2.0 * B_.zz() * lambda.value() - fr ) * t/lambda.value() );
  }

//Compute the updated A

volTensorField A = eigVecs_ & eigVals_ & eigVecs_.T() ;
