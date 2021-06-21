volTensorField gradU = fvc::grad(U) ;
//volSymmTensorField S = symm(gradU) ;
volTensorField S= gradU/2 + gradU.T()/2 ;

volTensorField R= gradU/2 - gradU.T()/2 ;

volTensorField B= ( eigVecs_.T() & S & eigVecs_ );

volTensorField D= ( eigVecs_.T() & R & eigVecs_ );


dimensionedTensor Itensor
( 
	"Identity", 
	dimensionSet(0, 0, 0, 0, 0, 0, 0), 
	tensor::I 
);

/*
volTensorField inter_M = ( eigVals_ & ( B + D) ) + ( (B - D) & eigVals_ );


forAll(eigVals_, cellI) //H é antisimétrica (demonstração no COllins) //H foi inicializado com a diagonal nula
 {
	tensor& H_ = H[cellI];
	tensor& inter_M_ = inter_M[cellI];
	tensor& vals = eigVals_[cellI];
	
	if ( abs(vals.yy() - vals.xx()) < 1e-8 ||  abs(vals.zz() - vals.xx()) < 1e-8 ||  abs(vals.zz() - vals.yy()) < 1e-8) {

  		H_.xy() = inter_M_.xy() / ( vals.yy() - vals.xx() + 1e-10);
		H_.xz() = inter_M_.xz() / ( vals.zz() - vals.xx() + 1e-10);
		H_.yz() = inter_M_.yz() / ( vals.zz() - vals.yy() + 1e-10);
	} 
	else {
		H_.xy() = inter_M_.xy() / ( vals.yy() - vals.xx() );
		H_.xz() = inter_M_.xz() / ( vals.zz() - vals.xx() );
		H_.yz() = inter_M_.yz() / ( vals.zz() - vals.yy() );
	}


	H_.yx() = - H_.xy();
	H_.zx() = - H_.xz();
	H_.zy() = - H_.yz();

	H_.xx() = 0;
	H_.yy() = 0;
	H_.zz() = 0;
 }
*/


forAll(eigVals_, cellI) //H é antisimétrica (demonstração no COllins) //H foi inicializado com a diagonal nula
 {
	tensor& H_ = H[cellI];
	tensor& vals = eigVals_[cellI];
	tensor& B_ = B[cellI];
	tensor& D_ = D[cellI];
	
	if ( abs(vals.yy() - vals.xx()) < 1e-8 ||  abs(vals.zz() - vals.xx()) < 1e-8 ||  abs(vals.zz() - vals.yy()) < 1e-8) {

  		H_.xy() = ( vals.xy() * ( B_.xy() + D_.xy()) ) + ( (B_.xy() - D_.xy()) * vals.xy() ) / ( vals.yy() - vals.xx() + 1e-10 );
		H_.xz() = ( vals.xz() * ( B_.xz() + D_.xz()) ) + ( (B_.xz() - D_.xz()) * vals.xz() ) / ( vals.zz() - vals.xx() + 1e-10 );
		H_.yz() = ( vals.yz() * ( B_.yz() + D_.yz()) ) + ( (B_.yz() - D_.yz()) * vals.yz() ) / ( vals.zz() - vals.yy() + 1e-10 );
	} 
	else {
		H_.xy() = ( vals.xy() * ( B_.xy() + D_.xy()) ) + ( (B_.xy() - D_.xy()) * vals.xy() ) / ( vals.yy() - vals.xx() );
		H_.xz() = ( vals.xz() * ( B_.xz() + D_.xz()) ) + ( (B_.xz() - D_.xz()) * vals.xz() ) / ( vals.zz() - vals.xx() );
		H_.yz() = ( vals.yz() * ( B_.yz() + D_.yz()) ) + ( (B_.yz() - D_.yz()) * vals.yz() ) / ( vals.zz() - vals.yy() );
	}


	H_.yx() = - H_.xy();
	H_.zx() = - H_.xz();
	H_.zy() = - H_.yz();

	H_.xx() = 0.;
	H_.yy() = 0.;
	H_.zz() = 0.;
 }

scalar deltat = runTime.deltaTValue();

 //Computing the multiplication directly using inv() funtion and matrix multiplication operator

eigVecs_ = Itensor & (inv( Itensor - deltat/2.0*somaHs* lambda/lambda.value()) & (Itensor + deltat/2.0 * somaHs * lambda/lambda.value())); //eigVecs0 = Itensor 


/*
forAll(eigVecs_, cellI) //Using the analytical solution given by collins
 {
	tensor& prod_ = prod[cellI];

	const float a = somaHs[cellI].zy() * deltat/2.0;
	const float b = somaHs[cellI].xz() * deltat/2.0;
	const float c = somaHs[cellI].yx() * deltat/2.0;

	const float a_squared = pow( a , 2.0);
	const float b_squared = pow( b , 2.0);
	const float c_squared = pow( c , 2.0);
 
	
	const float w = 1.0 / (1.0 + a_squared + b_squared + c_squared);

	prod_.xx() = 1 + a_squared - b_squared - c_squared;
	prod_.yy() = 1 - a_squared + b_squared - c_squared;
	prod_.zz() = 1 - a_squared - b_squared + c_squared;
	
	prod_.xy() = 2.0 * (a*b + c);
	prod_.yx() = 2.0 * (a*b - c);
	
	prod_.xz() = 2.0 * (a*c - b);
	prod_.zx() = 2.0 * (a*c + b);

	prod_.yz() = 2.0 * (b*c + a);
	prod_.zy() = 2.0 * (b*c - a);
	
	eigVecs_[cellI] =  tensor::I & ( w * prod_ ) ;

}
*/
somaHs +=  H - ( eigVecs_.T() & fvc::div(phi,eigVecs_) );

volScalarField Bxx = B.component(tensor::XX);
volScalarField Byy = B.component(tensor::YY);
volScalarField Bzz = B.component(tensor::ZZ);


/* // Numerically solving the differential equation for the temporal evolution of the principal deformations - using a kernel transformation with k()=ln():

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
*/

//Update of eigVals_

const float alpha = deltat*(L_squared.value()-3.0)/lambda.value();

volScalarField eigVals1 = eigVals_.component(tensor::XX);
volScalarField eigVals2 = eigVals_.component(tensor::YY);
volScalarField eigVals3 = eigVals_.component(tensor::ZZ);

volVectorField gradR_squared = fvc::grad(r_squared);

volScalarField U_gradR_squared = U & gradR_squared ; 

volScalarField gamma = eigVals1 * Foam::exp( Bxx * 2. * deltat * lambda/lambda.value() ) + eigVals2 * Foam::exp( Byy * 2. * deltat * lambda/lambda.value() ) + eigVals3 * Foam::exp( Bzz * 2. * deltat * lambda/lambda.value() ) - deltat * lambda/lambda.value() * U_gradR_squared + 3.*deltat/lambda.value();


//r²[n+1] must be computed first as a function of alpha and gamma
 
r_squared = 0.5 * ( gamma + alpha + L_squared - sqrt( pow( gamma + alpha + L_squared , 2.0 ) - 4.0 * gamma * L_squared) );

//f(r²[n+1]) must be now computed

volScalarField fr = ( L_squared - 3.0 ) / ( L_squared - r_squared );

//finally defining eig[i]:

volVectorField grad_eig1 = fvc::grad( eigVals_.component(tensor::XX) ) ;
volVectorField grad_eig2 = fvc::grad( eigVals_.component(tensor::YY) ) ;
volVectorField grad_eig3 = fvc::grad( eigVals_.component(tensor::ZZ) ) ;

volScalarField U_gradeig1 = U & grad_eig1 ;
volScalarField U_gradeig2 = U & grad_eig2 ;
volScalarField U_gradeig3 = U & grad_eig3 ;


//volScalarField eig1 = ( eig1 * Foam::exp(2.0 * deltat * Bxx * lambda/lambda.value() ) - deltat *U_gradeig1*lambda.value()/lambda + deltat/lambda.value() ) / ( 1 + deltat * fr / lambda.value() );

volScalarField eig1 = ( eigVals_.component(tensor::XX) * Foam::exp(2.0 * deltat * Bxx * lambda/lambda.value() ) - deltat *U_gradeig1*lambda/lambda.value() + deltat/lambda.value() ) / ( 1 + deltat * fr / lambda.value() );


volScalarField eig2 = ( eigVals_.component(tensor::YY) * Foam::exp(2.0 * deltat * Byy * lambda/lambda.value() ) - deltat *U_gradeig2*lambda/lambda.value() + deltat/lambda.value() ) / ( 1 + deltat * fr / lambda.value() );

volScalarField eig3 = ( eigVals_.component(tensor::ZZ) * Foam::exp(2.0 * deltat * Bzz * lambda/lambda.value() ) - deltat *U_gradeig3*lambda/lambda.value() + deltat/lambda.value() ) / ( 1 + deltat * fr / lambda.value() );

//Assing the values to the principal deformation matrix


//eigVals_.component(tensor::XX) = eig1; //the values can't be assigned in this manner, a loop is necessary:

forAll(eigVals_,cellI)  {

	if (eig1[cellI]<0) {   // setting negative eigenvalues to 0 to avoid Hadamard instabilities
	eigVals_[cellI].xx() = 0; }	
	else{
	eigVals_[cellI].xx() = eig1[cellI];	
	}
	
	if (eig2[cellI]<0) {
	eigVals_[cellI].yy() = 0; }	
	else{
	eigVals_[cellI].yy() = eig2[cellI];	
	}

	if (eig3[cellI]<0) {
	eigVals_[cellI].zz() = 0; }	
	else{
	eigVals_[cellI].zz() = eig3[cellI];	
	}

}

r_outro = eig1 + eig2 + eig3;

//Compute the conformation matrix

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );

// Compute the polymeric stress tensor:

tau_p = (nu_p/lambda) * symm( fr * A - Itensor);


tau_p.correctBoundaryConditions();
