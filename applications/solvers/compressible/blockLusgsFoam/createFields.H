#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pthermo(psiThermo::New(mesh));
psiThermo& thermo = pthermo();

autoPtr<gasProperties> pGasProps(gasProperties::New(thermo));
gasProperties& gasProps = pGasProps();

volScalarField& p = thermo.p();
volScalarField& T = thermo.T();
volScalarField& e = thermo.he();

//volScalarField h = e + p/thermo.rho();

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

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

/*volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(e + 0.5*magSqr(U))
);*/

volScalarField a
(
    IOobject
    (
        "a",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("a", dimVelocity, 0.0)
);

forAll(a, cell)
{
    a[cell] = gasProps.c(p[cell], T[cell]);
}  

convectiveFlux flux
(
    p,
    U,
    T,
    thermo
);

surfaceScalarField amaxSf
(
    "amaxSf",
    mag(fvc::interpolate(U) & mesh.Sf()) + mesh.magSf()*fvc::interpolate(a)
);

const surfaceScalarField& phi = flux.rhoFlux();
