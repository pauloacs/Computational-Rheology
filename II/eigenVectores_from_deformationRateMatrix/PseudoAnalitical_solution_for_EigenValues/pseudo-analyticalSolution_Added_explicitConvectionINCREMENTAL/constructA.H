volTensorField gradU = fvc::grad(U) ;
volSymmTensorField D = symm(gradU) ;



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
	//eigVals_[cellI].xx()=0;
   	//eigVals_[cellI].yy()=0;
   	//eigVals_[cellI].zz()=0;
	eigVals_[cellI].xx()=eival(0);
   	eigVals_[cellI].yy()=eival(1);
   	eigVals_[cellI].zz()=eival(2);
   }

 }


volScalarField f_r = (L2_ -3)/(L2_ - tr(eigVecs_ & eigVals_ & eigVecs_.T()) ) ;

volTensorField B= ( eigVecs_.T() & D & eigVecs_ );

//Evolute vals in time with analalytical solution for the transport of eigVals


scalar deltat = runTime.deltaTValue();


volScalarField eigx_conv = fvc::div( phi , eigVals_.component(tensor::XX) );
volScalarField eigy_conv = fvc::div( phi , eigVals_.component(tensor::YY) );
volScalarField eigz_conv = fvc::div( phi , eigVals_.component(tensor::ZZ) );


forAll(B, cellI)
 {
	tensor& B_ = B[cellI];
	scalar& fr = f_r[cellI];  	

	scalar termox = fr - 2.0 * B_.xx() * lambda.value();
	scalar termoy = fr - 2.0 * B_.yy() * lambda.value();
	scalar termoz = fr - 2.0 * B_.zz() * lambda.value();


	eigVals_[cellI].xx() *= Foam::exp( -termox * deltat/lambda.value() );
	eigVals_[cellI].xx() += (1.0 - lambda.value() * eigx_conv[cellI])/termox * (1.0 - Foam::exp( -termox * deltat/lambda.value() ) );

	eigVals_[cellI].yy() *= Foam::exp( -termoy * deltat/lambda.value() );
	eigVals_[cellI].yy() += (1.0 - lambda.value() * eigy_conv[cellI])/termoy * (1.0 - Foam::exp( -termoy * deltat/lambda.value() ) );
	
	eigVals_[cellI].zz() *= Foam::exp( -termoz * deltat/lambda.value() );
	eigVals_[cellI].zz() += (1.0 - lambda.value() * eigz_conv[cellI])/termoz * (1.0 - Foam::exp( -termoz * deltat/lambda.value() ) );
  }

//Compute the updated A

//volTensorField A = eigVecs_ & eigVals_ & eigVecs_.T() ;

f_r = (L2_-3)/(L2_ - tr(eigVecs_ & eigVals_ & eigVecs_.T()) ) ;


tau_p = (nu_p/lambda) * symm( f_r*(eigVecs_ & eigVals_ & eigVecs_.T()) - tensor::I);

tau_p.correctBoundaryConditions();
