volTensorField gradU = fvc::grad(U) ;

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


deltat = runTime.deltaT();


volScalarField Bxx = B.component(tensor::XX);
volScalarField Byy = B.component(tensor::YY);
volScalarField Bzz = B.component(tensor::ZZ);


//Update of eigVals_

dimensionedScalar alpha = deltat*(L_squared.value()-3.0)/lambda;

volScalarField eigVals1 = eigVals_.component(tensor::XX);
volScalarField eigVals2 = eigVals_.component(tensor::YY);
volScalarField eigVals3 = eigVals_.component(tensor::ZZ);

volVectorField gradR_squared = fvc::grad(r_squared);

volScalarField U_gradR_squared = U & gradR_squared ; 

volScalarField gamma = eigVals1 * Foam::exp( Bxx * 2. * deltat ) + eigVals2 * Foam::exp( Byy * 2. * deltat ) + eigVals3 * Foam::exp( Bzz * 2. * deltat ) - deltat * U_gradR_squared + 3.*deltat/lambda;


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


volScalarField eig1 = ( eigVals_.component(tensor::XX) * Foam::exp(2.0 * deltat * Bxx ) - deltat *U_gradeig1 + deltat/lambda ) / ( 1 + deltat * fr / lambda );


volScalarField eig2 = ( eigVals_.component(tensor::YY) * Foam::exp(2.0 * deltat * Byy ) - deltat *U_gradeig2 + deltat/lambda ) / ( 1 + deltat * fr / lambda );


volScalarField eig3 = ( eigVals_.component(tensor::ZZ) * Foam::exp(2.0 * deltat * Bzz ) - deltat *U_gradeig3 + deltat/lambda ) / ( 1 + deltat * fr / lambda );

//Assigning the values to the principal deformation matrix


forAll(eigVals_,cellI)  {

	if (eig1[cellI] < 0) {   // setting negative eigenvalues to 0 to avoid Hadamard instabilities
	eigVals_[cellI].xx() = 0; }	
	else {
	eigVals_[cellI].xx() = eig1[cellI];	
	}
	
	if (eig2[cellI] < 0) {
	eigVals_[cellI].yy() = 0; }	
	else {
	eigVals_[cellI].yy() = eig2[cellI];	
	}

	if (eig3[cellI] < 0) {
	eigVals_[cellI].zz() = 0; }	
	else {
	eigVals_[cellI].zz() = 1;	
	}

}


volTensorField inter_M = ( eigVals_ & ( B + D) ) + ( (B - D) & eigVals_ );


forAll(eigVals_, cellI) //H é antisimétrica (demonstração no COllins) //H foi inicializado com a diagonal nula
 {
	tensor& H_ = H[cellI];
	tensor& inter_M_ = inter_M[cellI];
	tensor& vals = eigVals_[cellI];

	H_.xy() = inter_M_.xy() / ( vals.yy() - vals.xx() + 1e-16);
	H_.xz() = inter_M_.xz() / ( vals.zz() - vals.xx() + 1e-16);
	H_.yz() = inter_M_.yz() / ( vals.zz() - vals.yy() + 1e-16);

	H_.yx() = - H_.xy();
	H_.zx() = - H_.xz();
	H_.zy() = - H_.yz();

	H_.xx() = 0;
	H_.yy() = 0;
	H_.zz() = 0;
 }


volTensorField term2 = Itensor + deltat/2.0 * somaHs;
volTensorField term1 = inv ( Itensor - deltat/2.0 * somaHs);

eigVecs_ = Itensor & term1 & term2;


volTensorField H_= H - ( eigVecs_.T() & fvc::div(phi,eigVecs_) );

somaHs +=  H_ ;


//Compute the conformation matrix

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );

// Compute the polymeric stress tensor:

tau_p = (nu_p/lambda) * symm( fr * A - Itensor);


tau_p.correctBoundaryConditions();
