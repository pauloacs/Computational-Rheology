Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties.lookup("nu")
);

dimensionedScalar nu_p
(
    "nu_p",
    transportProperties.lookup("nu_p")
);

dimensionedScalar lambda
(
    "lambda",
    transportProperties.lookup("lambda")
);

dimensionedScalar nu_stab
(
    "nu_stab",
    dimViscosity,
    transportProperties.lookup("nu_stab")
);

dimensionedScalar nu_stab_eigEqn
(
    "nu_stab_eigEqn",
    dimViscosity,
    transportProperties.lookup("nu_stab_eigEqn")
);

dimensionedScalar L2_
(
    "L2_",
    transportProperties.lookup("L2_")
);



Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volSymmTensorField tau_p
(
    IOobject
    (
        "tau_p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh

);

volTensorField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("H", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)

);

volTensorField eigVecs_
    (
        IOobject
        (
            "eigVecs_",
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
         extrapolatedCalculatedFvPatchField<tensor>::typeName
);

volTensorField eigVals_
    (
        IOobject
        (
            "eigVals_",
            U.time().timeName(),
            U.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedTensor
        (
                "I",
                dimless,
                pTraits<tensor>::I
        ),
        extrapolatedCalculatedFvPatchField<tensor>::typeName
);

volScalarField eig1
    (
        IOobject
        (
            "eig1",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
//        dimensionedScalar("eig1", dimensionSet(0,0,0,0,0,0,0), 0)

);

volScalarField eig2
    (
        IOobject
        (
            "eig2",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("eig2", dimensionSet(0,0,0,0,0,0,0), 0)
);
volScalarField eig3
    (
        IOobject
        (
            "eig3",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("eig3", dimensionSet(0,0,0,0,0,0,0), 0)
);

volScalarField fr
    (
        IOobject
        (
            "fr",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("fr", dimensionSet(0,0,0,0,0,0,0), 1)
);

volTensorField eigVecs0
    (
        IOobject
        (
            "eigVecs0",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("eigVecs0", dimensionSet(0,0,0,0,0,0,0), tensor::I)
);

volTensorField somaHs
    (
        IOobject
        (
            "somaHs",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
    	dimensionedTensor("somaHs", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
);


#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
