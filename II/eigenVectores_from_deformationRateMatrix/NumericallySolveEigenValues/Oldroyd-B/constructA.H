volTensorField gradU = fvc::grad(U) ;
volSymmTensorField D = symm(gradU) ; // corresponds to: volTensorField D= gradU/2 + gradU.T()/2 ;


const int t = runTime.value() * 1e4;
const int deltat = runTime.deltaTValue()*1e4;

const int deltat10 = 1 *deltat;

if ( t % deltat10  == 0.0 ) {
	// Eigen decomposition using a QR algorithm of Eigen library 
	Eigen::Matrix3d D_eig(Eigen::Matrix3d::Zero(3,3));
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;

	forAll(D, cellI)
	 {
	   if (tau_p[cellI].xx() < 100 ) {

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
		   //Eigen::Vector3d eival = eigSol.eigenvalues();
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
	  }
	 }
}

dimensionedTensor Itensor
( 
	"Identity", 
	dimensionSet(0, 0, 0, 0, 0, 0, 0), 
	tensor::I 
);


volTensorField B= ( eigVecs_.T() & D & eigVecs_ );

volScalarField Bxx = B.component(tensor::XX);
volScalarField Byy = B.component(tensor::YY);
volScalarField Bzz = B.component(tensor::ZZ);


fvScalarMatrix eigvalues1Eqn
            (
               fvm::ddt(eig1) + fvm::div(phi,eig1)
             ==
               2.0 *Bxx  +  (  Foam::exp(-eig1) - 1.0  )/lambda 
            );

fvScalarMatrix eigvalues2Eqn
            (
                fvm::ddt(eig2) + fvm::div(phi,eig2)
             ==
               2.0 *Byy  +  ( Foam::exp(-eig2) - 1.0 )/lambda 
            );

fvScalarMatrix eigvalues3Eqn
            (
                fvm::ddt(eig3) + fvm::div(phi,eig3)
             ==
               2.0 *Bzz  +  ( Foam::exp(-eig3) - 1.0 )/lambda
            );


eigvalues1Eqn.relax();
eigvalues1Eqn.solve();

eigvalues2Eqn.relax();
eigvalues2Eqn.solve();

eigvalues3Eqn.relax();
eigvalues3Eqn.solve();


forAll(eigVals_, cellI)
 {
	eigVals_[cellI].xx() = Foam::exp( eig1[cellI] );
	eigVals_[cellI].yy() = Foam::exp( eig2[cellI] );
	eigVals_[cellI].zz() = Foam::exp( eig3[cellI] );
 }


volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );
/*
forAll(A, cellI)
 {
        A[cellI].zz() = 0.0 ;
 }
*/
//A.write();


tau_p = (nu_p/lambda) * symm( A - Itensor);
/*
forAll(tau_p, cellI)
 {
        tau_p[cellI].zz() = 0.0 ;
 }
*/
tau_p.correctBoundaryConditions();

//--------------------------------------------------------------------------------------//

