#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

TwoFluidFoam::twoFluid fluidSystem(mesh);

twoFluidConservative& conservative = fluidSystem.conservative();

twoFluidConvectiveFlux twoFluidFlux
(
    fluidSystem.p(),
    fluidSystem.alpha(),
    fluidSystem.U1(),
    fluidSystem.U2(),
    fluidSystem.T1(),
    fluidSystem.T2(),
    fluidSystem.thermo1(),
    fluidSystem.thermo2()
);


volScalarField& p = fluidSystem.p();
const volScalarField& alpha1 = fluidSystem.alpha1();
volScalarField& alpha2 = fluidSystem.alpha2();
volVectorField& U1 = fluidSystem.U1();
volVectorField& U2 = fluidSystem.U2();
volScalarField& T1 = fluidSystem.T1();
volScalarField& T2 = fluidSystem.T2();

rhoThermo& thermo1 = fluidSystem.thermo1();
rhoThermo& thermo2 = fluidSystem.thermo2();



surfaceScalarField amaxSf
(
    "amaxSf",
    fluidSystem.amaxSf()
);

surfaceScalarField alphaRhoPhi2
(
    IOobject
    (
        "alphaRhoPhi2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(thermo2.rho())*fvc::interpolate(alpha2)*fvc::flux(U2)
);

autoPtr<TwoFluidFoam::dragModel> pDrag = TwoFluidFoam::dragModel::New(fluidSystem);
TwoFluidFoam::dragModel& drag = pDrag();

autoPtr<liquidProperties> pLiquidProps = liquidProperties::New(fluidSystem.thermo2().lookupOrDefault<word>("liquidProperties", "H2O"));
liquidProperties& liquidProps = pLiquidProps();
autoPtr<saturationCurve> pSaturation = saturationCurve::New(fluidSystem);
saturationCurve& saturation = pSaturation();


autoPtr<Foam::WetSteam::condensationModel> pCondensation = Foam::WetSteam::condensationModel::New
(
    fluidSystem.alpha2(),
    fluidSystem.thermo2().rho(),
    fluidSystem.U2(),
    alphaRhoPhi2,
    fluidSystem.thermo1(),
    fluidSystem.thermo2(),
    saturation,
    liquidProps
);
Foam::WetSteam::condensationModel& condensation = pCondensation();


//Pro ukladani
volScalarField rho1
(
    IOobject
    (
        "rho1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo1.rho()
);

volScalarField rho2
(
    IOobject
    (
        "rho2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo2.rho()
);

volScalarField e2
(
    IOobject
    (
        "e2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo2.he()
);

