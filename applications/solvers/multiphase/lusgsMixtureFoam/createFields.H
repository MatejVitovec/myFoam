//#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

compressibleMixture fluid(mesh);

compressibleMixtureConvectiveFlux mixtureFlux
(
    fluid.p(),
    fluid.U(),
    fluid.T(),
    fluid.w(),
    fluid.thermo1(),
    fluid.thermo2()
);

volScalarField& p = fluid.p();
volVectorField& U = fluid.U();
volScalarField& T = fluid.T();
volScalarField& w = fluid.w();

volScalarField alpha = fluid.alpha();

rhoThermo& thermo1 = fluid.thermo1();
rhoThermo& thermo2 = fluid.thermo2();

const volScalarField& a = fluid.a();

const volScalarField& ev = thermo1.he();
const volScalarField& el = thermo2.he();

const volScalarField& rhov = thermo1.rho();
const volScalarField& rhol = thermo2.rho();

volScalarField rho
(
    IOobject
    (
        "rho",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fluid.rho()
);

surfaceScalarField amaxSf
(
    "amaxSf",
    fluid.amaxSf()
);


volScalarField dp
(
    IOobject
    (
        "dp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dp", dimPressure, 0.0)
);

volVectorField dU
(
    IOobject
    (
        "dU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("dU", dimVelocity, vector(0, 0, 0))
);

volScalarField dT
(
    IOobject
    (
        "dT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dT", dimTemperature, 0.0)
);

volScalarField dw
(
    IOobject
    (
        "dw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dw", dimless, 0.0)
);
