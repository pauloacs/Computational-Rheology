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



fvScalarMatrix eigvalues1Eqn
            (
               fvm::ddt(eig1) + fvm::div(phi,eig1)
             ==
               2.0 *Bxx  +  (  Foam::exp(-eig1) - 1  )/lambda 
            );

fvScalarMatrix eigvalues2Eqn
            (
                fvm::ddt(eig2) + fvm::div(phi,eig2)
             ==
               2.0 *Byy  +  ( Foam::exp(-eig2) - 1 )/lambda 
            );

fvScalarMatrix eigvalues3Eqn
            (
                fvm::ddt(eig3) + fvm::div(phi,eig3)
             ==
               2.0 *Bzz  +  ( Foam::exp(-eig3) - 1 )/lambda
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


volScalarField H21 = H.component(tensor::YX);
volScalarField H12 = - H21;

volScalarField H13 = H.component(tensor::XZ);
volScalarField H31 = - H13;


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



fvScalarMatrix eigvecs21Eqn
            (
               fvm::ddt(v21)  - fvm::laplacian(nu_stab1,v21) //+ fvm::div(phi,v21)
             == H21 * v22 - u_nablaU.component(tensor::YX)  - nu_stab1 * fvc::div(fvc::grad(v21))
                
            );



eigvecs21Eqn.relax();
eigvecs21Eqn.solve();

forAll(eigVecs_, cellI) {

	if (v21[cellI] >1 ) {
	
		v21[cellI] = 1.0;
		v11[cellI] = 0.0;
	}
	else if ( v21[cellI] < -1 ){

		v21[cellI] = -1.0;
		v11[cellI] = 0.0;
	}
}



fvScalarMatrix eigvecs22Eqn
            (
               
               fvm::ddt(v22)  - fvm::laplacian(nu_stab1,v22) //+ fvm::div(phi,v23)
             == H12 * v21 - u_nablaU.component(tensor::YY)   - nu_stab1 * fvc::div(fvc::grad(v22))
               
            );


eigvecs22Eqn.relax();
eigvecs22Eqn.solve();


forAll(eigVecs_, cellI) {

	if (v22[cellI] >1 ) {
	
		v22[cellI] = 1.0;
		v12[cellI] = 0.0;
	}
	else if ( v22[cellI] < -1 ){

		v22[cellI] = -1.0;
		v12[cellI] = 0.0;
	}

}

fvScalarMatrix eigvecs11Eqn
            (
               fvm::ddt(v11)  - fvm::laplacian(nu_stab1,v11) //+ fvm::div(phi,v11)
             ==  H21 * v12 - u_nablaU.component(tensor::XX) - nu_stab1 * fvc::div(fvc::grad(v11)) 
                
            );

eigvecs11Eqn.relax();
eigvecs11Eqn.solve();


forAll(eigVecs_, cellI) {

	if (v11[cellI] >1 ) {
	
		v11[cellI] = 1.0;
		v21[cellI] = 0.0;
	}
	else if ( v11[cellI] < -1 ){

		v11[cellI] = -1.0;
		v21[cellI] = 0.0;
	}

}



fvScalarMatrix eigvecs12Eqn
            (
               fvm::Sp(v11,v12) //- fvm::laplacian(nu_stab1,v13)
             == - v21 *v22 //- nu_stab1 * fvc::div(fvc::grad(v13))
            );

eigvecs12Eqn.relax();
eigvecs12Eqn.solve();


forAll(eigVecs_, cellI) {

	if (v12[cellI] >1 ) {
	
		v12[cellI] = 1.0;
		v22[cellI] = 0.0;
	}
	else if ( v12[cellI] < -1 ){

		v12[cellI] = -1.0;
		v22[cellI] = 0.0;
	}

}
	

//CONSTRUCT U

forAll(eigVecs_, cellI) {
	
		tensor& eigVecs = eigVecs_[cellI];

		eigVecs.xx() = v11[cellI];
		eigVecs.yx() = v21[cellI];
		eigVecs.zx() = 0;

		eigVecs.xy() = v12[cellI];
		eigVecs.yy() = v22[cellI];
		eigVecs.zy() = 0;
	
		eigVecs.xz() = 0;
		eigVecs.yz() = 0;
		eigVecs.zz() = 1.0;
	
}

volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );

tau_p = (nu_p/lambda) * symm( A - Itensor);

tau_p.correctBoundaryConditions();

