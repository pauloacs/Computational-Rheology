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

deltat = runTime.deltaT();

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


volTensorField inter_M = (eigVals_ & ( B + D)) + ((B - D) & eigVals_);


forAll(eigVals_, cellI) //H is skew symmetric
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



dimensionedTensor U0
( 
	"U0", 
	dimensionSet(0, 0, 0, 0, 0, 0, 0), 
	tensor(1,0,0,0,0,1,0,1,0)
);

volTensorField term2 = Itensor + deltat/2.0 * somaHs;
volTensorField term1 = inv ( Itensor - deltat/2.0 * somaHs);


eigVecs_ = Itensor & term1 & term2;

forAll(eigVecs_, cellI) 
 {
	if (eigVecs_[cellI].xx() > 1.0 ) {
		
		eigVecs_[cellI].xx() = 1.0;
		eigVecs_[cellI].yx() = 0;

	}
	else if (eigVecs_[cellI].xx() < -1.0 ) {
		
		eigVecs_[cellI].xx() = - 1.0;
		eigVecs_[cellI].yx() = 0;

	}

	if (eigVecs_[cellI].yx() > 1.0 ) {
		
		eigVecs_[cellI].yx() = 1.0;
		eigVecs_[cellI].xx() = 0;

	}
	else if (eigVecs_[cellI].yx() < -1.0 ) {
		
		eigVecs_[cellI].yx() = - 1.0;
		eigVecs_[cellI].xx() = 0;

	}

	if (eigVecs_[cellI].xy() > 1.0 ) {
		
		eigVecs_[cellI].xy() = 1.0;
		eigVecs_[cellI].yy() = 0;

	} 
	else if (eigVecs_[cellI].xy() < -1.0 ) {
		
		eigVecs_[cellI].xy() = - 1.0;
		eigVecs_[cellI].yy() = 0;

	}

	if (eigVecs_[cellI].yy() > 1.0 ) {
		
		eigVecs_[cellI].yy() = 1.0;
		eigVecs_[cellI].xy() = 0;

	}
	else if (eigVecs_[cellI].yy() < -1.0 ) {
		
		eigVecs_[cellI].yy() = - 1.0;
		eigVecs_[cellI].xy() = 0;

	}

 }

//DOING U & GRAD(EIGVECS_) instead of fvc::div(phi,eigVecs_)

volVectorField grad_eigVec11 = fvc::grad( eigVecs_.component(tensor::XX) );
volVectorField grad_eigVec21 = fvc::grad( eigVecs_.component(tensor::YX) );
volVectorField grad_eigVec31 = fvc::grad( eigVecs_.component(tensor::ZX) );

volVectorField grad_eigVec12 = fvc::grad( eigVecs_.component(tensor::XY) );
volVectorField grad_eigVec22 = fvc::grad( eigVecs_.component(tensor::YY) );
volVectorField grad_eigVec32 = fvc::grad( eigVecs_.component(tensor::ZY) );

volVectorField grad_eigVec13 = fvc::grad( eigVecs_.component(tensor::XZ) );
volVectorField grad_eigVec23 = fvc::grad( eigVecs_.component(tensor::YZ) );
volVectorField grad_eigVec33 = fvc::grad( eigVecs_.component(tensor::ZZ) );


forAll(eigVecs_, cellI) 
 {
	tensor& u_nablaU_ = u_nablaU[cellI];

	vector& U_ = U[cellI];
	vector v1(0,0,0);
	vector v2(0,0,0);
	vector v3(0,0,0);
	vector v4(0,0,0);
	vector v5(0,0,0);
	vector v6(0,0,0);
	vector v7(0,0,0);
	vector v8(0,0,0);
	vector v9(0,0,0);

	v1.x() = grad_eigVec11[cellI].x();
	v1.y() = grad_eigVec21[cellI].x();
	v1.z() = grad_eigVec31[cellI].x();
	
	v2.x() = grad_eigVec11[cellI].y();
	v2.y() = grad_eigVec21[cellI].y();
	v2.z() = grad_eigVec31[cellI].y();

	v3.x() = grad_eigVec11[cellI].z();
	v3.y() = grad_eigVec21[cellI].z();
	v3.z() = grad_eigVec31[cellI].z();

	v4.x() = grad_eigVec12[cellI].x();
	v4.y() = grad_eigVec22[cellI].x();
	v4.z() = grad_eigVec32[cellI].x();

	v5.x() = grad_eigVec12[cellI].y();
	v5.y() = grad_eigVec22[cellI].y();
	v5.z() = grad_eigVec32[cellI].y();

	v6.x() = grad_eigVec12[cellI].z();
	v6.y() = grad_eigVec22[cellI].z();
	v6.z() = grad_eigVec32[cellI].z();

	v7.x() = grad_eigVec13[cellI].x();
	v7.y() = grad_eigVec23[cellI].x();
	v7.z() = grad_eigVec33[cellI].x();

	v8.x() = grad_eigVec13[cellI].y();
	v8.y() = grad_eigVec23[cellI].y();
	v8.z() = grad_eigVec33[cellI].y();

	v9.x() = grad_eigVec13[cellI].z();
	v9.y() = grad_eigVec23[cellI].z();
	v9.z() = grad_eigVec33[cellI].z();

	u_nablaU_.xx() = v1 & U_;
	u_nablaU_.xy() = v2 & U_;
	u_nablaU_.xz() = v3 & U_;
	
	u_nablaU_.yx() = v4 & U_;
	u_nablaU_.yy() = v5 & U_;
	u_nablaU_.yz() = v6 & U_;

	u_nablaU_.zx() = v7 & U_;
	u_nablaU_.zy() = v8 & U_;
	u_nablaU_.zz() = v9 & U_;

 }

H_ast = H - ( eigVecs_.T() & u_nablaU);

somaHs += H_ast;

//OR H_ast = H - ( eigVecs_.T() & fvc::div(phi,eigVecs_) );
//

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );

tau_p = (nu_p/lambda) * symm( A - Itensor);

tau_p.correctBoundaryConditions();

