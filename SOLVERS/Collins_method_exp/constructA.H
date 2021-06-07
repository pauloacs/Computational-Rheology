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

volTensorField inter_M = (eigVals_& ( B + D)) + ((B - D) & eigVals_);

forAll(eigVals_, cellI) //H é antisimétrica (demonstração no COllins) //H foi inicializado com a diagonal nula
 {
	tensor& H_ = H[cellI];
	tensor& inter_M_ = inter_M[cellI];
	tensor& vals = eigVals_[cellI];
	
	H_.xy() = inter_M_.xy() / ( vals.yy() - vals.xx() + 1e8);
	H_.xz() = inter_M_.xz() / ( vals.zz() - vals.xx() + 1e8);
	H_.yz() = inter_M_.yz() / ( vals.zz() - vals.yy() + 1e8);

	H_.yx() = - H_.xy();
	H_.zx() = - H_.xz();
	H_.zy() = - H_.yz();

	H_.xx() = 0;
	H_.yy() = 0;
	H_.zz() = 0;
 }

scalar deltat = runTime.deltaTValue();

somaHs +=  H - ( eigVecs_.T() & fvc::div(phi,eigVecs_) );

eigVecs_ = Itensor & (inv( Itensor - deltat/2.0*somaHs* lambda/lambda.value()) & (Itensor + deltat/2.0*somaHs* lambda/lambda.value())); //eigVecs0 = Itensor 


forAll(fr, cellI) 
 {
	fr[cellI] = (L2_.value() - 3.0)/( L2_.value() - (eigVals_[cellI].xx() + eigVals_[cellI].zz() + eigVals_[cellI].yy() ) );
 }

volScalarField Bxx = B.component(tensor::XX);
volScalarField Byy = B.component(tensor::YY);
volScalarField Bzz = B.component(tensor::ZZ);


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


volSymmTensorField A= symm( eigVecs_ & eigVals_ & eigVecs_.T() );


forAll(fr, cellI) 
 {
	fr[cellI] = (L2_.value() - 3.0)/( L2_.value() - (eigVals_[cellI].xx() + eigVals_[cellI].zz() + eigVals_[cellI].yy() ) );
 }

tau_p = (nu_p/lambda) * symm( fr * A - Itensor);

tau_p.correctBoundaryConditions();
