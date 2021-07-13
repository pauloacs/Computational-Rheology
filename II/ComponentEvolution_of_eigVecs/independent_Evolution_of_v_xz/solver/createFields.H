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

dimensionedScalar nu_stab1
(
    "nu_stab1",
    transportProperties.lookup("nu_stab1")
);

dimensionedScalar L_squared
(
    "L_squared",
    transportProperties.lookup("L_squared")
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

volScalarField r_squared
    (
        IOobject
        (
            "r_squared",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("r_squared", dimensionSet(0,0,0,0,0,0,0), 3.0)
);

volTensorField prod
    (
        IOobject
        (
            "prod",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("prod", dimensionSet(0,0,0,0,0,0,0), tensor::I)
);

volTensorField H_
    (
        IOobject
        (
            "H_",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
    	dimensionedTensor("H_", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
);

volTensorField H__
    (
        IOobject
        (
            "H__",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
    	dimensionedTensor("H__", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
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

volTensorField u_nablaU
    (
        IOobject
        (
            "u_nablaU",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
    	dimensionedTensor("u_nablaU", dimensionSet(0,0,-1,0,0,0,0), tensor::zero)
);

volScalarField v11
    (
        IOobject
        (
            "v11",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
);
volScalarField v21
    (
        IOobject
        (
            "v21",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
);

volScalarField v13
    (
        IOobject
        (
            "v13",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
        //dimensionedScalar("v12", dimensionSet(0,0,0,0,0,0,0), 0.0)
);

volScalarField v23
    (
        IOobject
        (
            "v23",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
);



#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
