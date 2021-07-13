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

scalar deltat = runTime.deltaTValue();

volScalarField Bxx = B.component(tensor::XX);
volScalarField Byy = B.component(tensor::YY);
volScalarField Bzz = B.component(tensor::ZZ);


//---------------------------------<DETERMINAR fr>--------------------------------------------------//


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

//-------------------------------</DETERMINAR fr>------------------------------------------------//


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

somaHs +=  H_;


volTensorField inter_M = (eigVals_ & ( B + D)) + ((B - D) & eigVals_);

forAll(eigVals_, cellI) 
 {
	tensor& H_ = H[cellI];
	tensor& inter_M_ = inter_M[cellI];
	tensor& vals = eigVals_[cellI];

	H_.xy() = inter_M_.xy() / ( vals.yy() - vals.xx() + 1e-20 );
	H_.xz() = inter_M_.xz() / ( vals.zz() - vals.xx() + 1e-20 );
	H_.yz() = inter_M_.yz() / ( vals.zz() - vals.yy() + 1e-20 );	
	
	H_.yx() = - H_.xy();
	H_.zx() = - H_.xz();
	H_.zy() = - H_.yz();

	H_.xx() = 0;
	H_.yy() = 0;
	H_.zz() = 0;
 }

/* //relation for (I-A)⁻¹(I+A) in the article
forAll(eigVecs_, cellI) 
 {
	tensor& prod_ = prod[cellI];

	const float a = somaHs[cellI].zy() * deltat.value()/2.0;
	const float b = somaHs[cellI].xz() * deltat.value()/2.0;
	const float c = somaHs[cellI].yx() * deltat.value()/2.0;

	const float a_squared = pow( a , 2.0);
	const float b_squared = pow( b , 2.0);
	const float c_squared = pow( c , 2.0);
 
	
	const float w = (1.0 + a_squared + b_squared + c_squared);

	prod_.xx() = 1 + a_squared - b_squared - c_squared;
	prod_.yy() = 1 - a_squared + b_squared - c_squared;
	prod_.zz() = 1 - a_squared - b_squared + c_squared;
	
	prod_.xy() = 2.0 * (a * b + c );
	prod_.yx() = 2.0 * (a * b - c );
	
	prod_.xz() = 2.0 * (a * c - b );
	prod_.zx() = 2.0 * (a * c + b );

	prod_.yz() = 2.0 * (b * c + a);
	prod_.zy() = 2.0 * (b * c - a );
	
	tensor U0(1,0,0,0,0,1,0,1,0);
	eigVecs_[cellI] =  U0 & ( prod_ / w ) ;

}
*/



eigVecs_ = ( inv( Itensor - deltat/2.0*somaHs* lambda/lambda.value() ) ) & ( Itensor + deltat/2.0 * somaHs * lambda/lambda.value() ); //eigVecs0 = Itensor

H_ = H - ( eigVecs_.T() & fvc::div(phi,eigVecs_) );


/*//using U & GRAD(EIGVECS_) 

volScalarField eigVecs11 = eigVecs_.component(tensor::XX);
volScalarField eigVecs22 = eigVecs_.component(tensor::YY);
volScalarField eigVecs33 = eigVecs_.component(tensor::ZZ);

volScalarField eigVecs12 = eigVecs_.component(tensor::XY);
volScalarField eigVecs13 = eigVecs_.component(tensor::XZ);

volScalarField eigVecs21 = eigVecs_.component(tensor::YX);
volScalarField eigVecs23 = eigVecs_.component(tensor::YZ);

volScalarField eigVecs31 = eigVecs_.component(tensor::ZX);
volScalarField eigVecs32 = eigVecs_.component(tensor::ZY);

volScalarField U_grad_eigVec11 = U & fvc::grad(eigVecs11);
volScalarField U_grad_eigVec22 = U & fvc::grad(eigVecs22);
volScalarField U_grad_eigVec33 = U & fvc::grad(eigVecs33);

volScalarField U_grad_eigVec12 = U & fvc::grad(eigVecs12);
volScalarField U_grad_eigVec13 = U & fvc::grad(eigVecs13);

volScalarField U_grad_eigVec21 = U & fvc::grad(eigVecs21);
volScalarField U_grad_eigVec23 = U & fvc::grad(eigVecs23);

volScalarField U_grad_eigVec31 = U & fvc::grad(eigVecs31);
volScalarField U_grad_eigVec32 = U & fvc::grad(eigVecs32);


forAll(U_grad_eigVec, cellI) 
 {
	tensor& U_grad_eigVec_ = U_grad_eigVec[cellI];

	U_grad_eigVec_.xx() = U_grad_eigVec11[cellI];
	U_grad_eigVec_.yy() = U_grad_eigVec22[cellI];
	U_grad_eigVec_.zz() = U_grad_eigVec33[cellI];

	U_grad_eigVec_.xy() = U_grad_eigVec12[cellI];
	U_grad_eigVec_.xz() = U_grad_eigVec13[cellI];

	U_grad_eigVec_.yx() = U_grad_eigVec21[cellI];
	U_grad_eigVec_.yz() = U_grad_eigVec23[cellI];

	U_grad_eigVec_.zx() = U_grad_eigVec31[cellI];
	U_grad_eigVec_.zy() = U_grad_eigVec32[cellI];
 }


//somaHs +=  H - U_grad_eigVec ;

//H__ = H - ( eigVecs_.T() & U_grad_eigVec );
*/

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );
//A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );


tau_p = (nu_p/lambda) * symm( fr * A - Itensor);

tau_p.correctBoundaryConditions();

