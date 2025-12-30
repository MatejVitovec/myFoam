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
volScalarField rho = fluid.rho();

rhoThermo& thermo1 = fluid.thermo1();
rhoThermo& thermo2 = fluid.thermo2();

const volScalarField& a = fluid.a();

const volScalarField& e1 = thermo1.he();
const volScalarField& e2 = thermo2.he();

const volScalarField& rho1 = thermo1.rho();
const volScalarField& rho2 = thermo2.rho();

surfaceScalarField amaxSf
(
    "amaxSf",
    fluid.amaxSf()
);

const volScalarField alphaOne
(
    IOobject
    (
        "alphaOne",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaOne", dimless, 1.0)
);

saturation satur = saturation(fluid.subDict("condensation"), thermo1);

autoPtr<Foam::WetSteam::condensationModel> pCondensation = Foam::WetSteam::condensationModel::New
(
    alpha,
    alphaOne,
    rho,
    fluid.U(),
    mixtureFlux.rhoFlux(),
    fluid.thermo1(),
    fluid.thermo2(),
    satur,
    fluid
);
Foam::WetSteam::condensationModel& condensation = pCondensation();


volScalarField condensationMassSource
(
    IOobject
    (
        "mDotCondensation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    condensation.condensationRateMassSource()
);

volScalarField nucleationMassSource
(
    IOobject
    (
        "mDotNucleation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    condensation.nucleationRateMassSource()
);

volScalarField growthMassSource
(
    IOobject
    (
        "mDotGrowth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    condensation.growthRateMassSource()
);

volScalarField s1
(
    IOobject
    (
        "s.1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("s", dimEnergy/dimMass/dimTemperature, 0.0)
);

volScalarField h1
(
    IOobject
    (
        "h.1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("h1", dimEnergy/dimMass, 0.0)
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
