volTensorField gradU = fvc::grad(U) ;
volSymmTensorField D = symm(gradU) ; // corresponds to: volTensorField D= gradU/2 + gradU.T()/2 ;



const int deltat = runTime.deltaTValue()*1e4;
const int t = runTime.value() * 1e4;
scalar count_cells = 0;


// Eigen decomposition using a QR algorithm of Eigen library 
Eigen::Matrix3d D_eig(Eigen::Matrix3d::Zero(3,3));
Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigSol;

forAll(D, cellI)
 {
    if ( abs ( mag( tau_p[cellI] ) ) < 30.0 ) {
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

const int deltat21 = 21.0 * deltat;

if ( t % deltat21  == 0.0 ) {
 
   count = 0;
   forAll(eigVecs_, cellI) {

	   if ( abs( eigVecs_[cellI].xx() - eigVecs_old[cellI].xx() )  > 0.4 ) {
		
		count+= 1;
		
	   }
    	      count_cells+=1;
   } 

   cout << "Cells with significative variation in eigen vectors : " << count <<"\n" ;
   cout << "Total number of cells" << count_cells <<"\n";
   eigVecs_old = eigVecs_;
}

cout << "got here \n";


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

//--------------------------------------------------------------------------------------//


const float alpha = deltat*(L_squared.value()-3.0)/lambda.value();

volScalarField eigVals1 = eigVals_.component(tensor::XX);
volScalarField eigVals2 = eigVals_.component(tensor::YY);
volScalarField eigVals3 = eigVals_.component(tensor::ZZ);

volVectorField gradR_squared = fvc::grad(r_squared);

volScalarField U_gradR_squared = U & gradR_squared ; 

volScalarField gamma = eigVals1 * Foam::exp( Bxx * 2. * deltat * lambda/lambda.value() ) + eigVals2 * Foam::exp( Byy * 2. * deltat * lambda/lambda.value() ) + eigVals3 * Foam::exp( Bzz * 2. * deltat * lambda/lambda.value() ) - deltat * lambda/lambda.value() * U_gradR_squared + 3.*deltat/lambda.value();


//r²[n+1] must be computed first as a function of alpha and gamma
 
r_squared = 0.5 * ( gamma + alpha + L_squared - sqrt( pow( gamma + alpha + L_squared , 2.0 ) - 4.0 * gamma * L_squared) );

//f(r[n+1]) must be now computed

volScalarField fr = ( L_squared - 3.0 ) / ( L_squared - r_squared );

//--------------------------------------------------------------------------------------//


fvScalarMatrix eigvalues1Eqn
            (
               fvm::ddt(eig1) + fvm::div(phi,eig1)
             ==
               2.0 *Bxx  +  (  Foam::exp(-eig1) - fr  )/lambda 
            );

fvScalarMatrix eigvalues2Eqn
            (
                fvm::ddt(eig2) + fvm::div(phi,eig2)
             ==
               2.0 *Byy  +  ( Foam::exp(-eig2) - fr )/lambda 
            );

fvScalarMatrix eigvalues3Eqn
            (
                fvm::ddt(eig3) + fvm::div(phi,eig3)
             ==
               2.0 *Bzz  +  ( Foam::exp(-eig3) - fr )/lambda
            );

if (  ( count > 0.01 * count_cells ) || ( t < 1000*deltat ) ) {

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

}

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );
/*
forAll(A, cellI)
 {
        A[cellI].zz() = 0.0 ;
 }
*/
//A.write();


tau_p = (nu_p/lambda) * symm( fr * A - Itensor);
/*
forAll(tau_p, cellI)
 {
        tau_p[cellI].zz() = 0.0 ;
 }
*/
tau_p.correctBoundaryConditions();
