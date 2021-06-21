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

forAll(eigVals_, cellI) //H é antisimétrica (demonstração no COllins) //H foi inicializado com a diagonal nula
 {
	tensor& H_ = H[cellI];
	tensor& vals = eigVals_[cellI];
	tensor& B_ = B[cellI];
	tensor& D_ = D[cellI];
	
	
	if ( abs(vals.yy() - vals.xx()) < 1e-25 ||  abs(vals.zz() - vals.xx()) < 1e-25 ||  abs(vals.zz() - vals.yy()) < 1e-25) {

  		H_.xy() = ( vals.xy() * ( B_.xy() + D_.xy()) ) + ( (B_.xy() - D_.xy()) * vals.xy() ) / ( vals.yy() - vals.xx() + 1e-25 );
		H_.xz() = ( vals.xz() * ( B_.xz() + D_.xz()) ) + ( (B_.xz() - D_.xz()) * vals.xz() ) / ( vals.zz() - vals.xx() + 1e-25 );
		H_.yz() = ( vals.yz() * ( B_.yz() + D_.yz()) ) + ( (B_.yz() - D_.yz()) * vals.yz() ) / ( vals.zz() - vals.yy() + 1e-25 );
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


//eigVecs_ = Itensor & (inv( Itensor - deltat/2.0*somaHs* lambda/lambda.value()) & (Itensor + deltat/2.0 * somaHs * lambda/lambda.value())); //eigVecs0 = Itensor 

forAll(eigVecs_, cellI) 
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

somaHs +=  H - ( eigVecs_.T() & fvc::div(phi,eigVecs_) );

/*
forAll(fr, cellI) 
 {
	fr[cellI] = (L2_.value() - 3.0)/( L2_.value() - (eigVals_[cellI].xx() + eigVals_[cellI].zz() + eigVals_[cellI].yy() ) );
 }
*/

volScalarField Bxx = B.component(tensor::XX);
volScalarField Byy = B.component(tensor::YY);
volScalarField Bzz = B.component(tensor::ZZ);


//--------------------------------------------------------------------------------------//

//Compute f(r[n+1])

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

/*
forAll(fr, cellI) 
 {
	fr[cellI] = (L2_.value() - 3.0)/( L2_.value() - (eigVals_[cellI].xx() + eigVals_[cellI].zz() + eigVals_[cellI].yy() ) );
 }
*/

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );

tau_p = (nu_p/lambda) * symm( fr * A - Itensor);

tau_p.correctBoundaryConditions();

