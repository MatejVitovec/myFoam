#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo1
(
    psiThermo::New(mesh)
);
psiThermo& thermo1 = pThermo1();

autoPtr<psiThermo> pThermo2
(
    psiThermo::New(mesh)
);
psiThermo& thermo2 = pThermo2();

volScalarField& p = thermo1.p();
volScalarField& T1 = thermo1.T();
volScalarField& T2 = thermo2.T(); //TODO

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField pInt
(
    IOobject
    (
        "pInt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

//Info<< "Reading field U1\n" << endl;
volVectorField U1
(
    IOobject
    (
        "U1",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField U2
(
    IOobject
    (
        "U2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// Conservative

volScalarField alphaRho1
(
    IOobject
    (
        "alphaRho1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha*thermo1.rho()
);

volScalarField alphaRho2
(
    IOobject
    (
        "alphaRho2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (1.0 - alpha)*thermo2.rho()
);

volVectorField alphaRhoU1
(
    IOobject
    (
        "alphaRhoU1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaRho1*U1
);

volVectorField alphaRhoU2
(
    IOobject
    (
        "alphaRhoU2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaRho2*U2
);

volScalarField alphaRhoE1
(
    IOobject
    (
        "alphaRhoE1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaRho1*((thermo1.he() + 0.5*magSqr(U1)) - p/thermo1.rho())
);

volScalarField alphaRhoE2
(
    IOobject
    (
        "alphaRhoE2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaRho2*((thermo2.he() + 0.5*magSqr(U2)) - p/thermo2.rho())
);


twoFluidConvectiveFlux dbFlux
(
    p,
    alpha,
    U1,
    U2,
    T1,
    T2,
    thermo1,
    thermo2
);
