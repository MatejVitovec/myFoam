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
volScalarField& alpha = fluidSystem.alpha();
volVectorField& U1 = fluidSystem.U1();
volVectorField& U2 = fluidSystem.U2();
volScalarField& T1 = fluidSystem.T1();
volScalarField& T2 = fluidSystem.T2();

rhoThermo& thermo1 = fluidSystem.thermo1();
rhoThermo& thermo2 = fluidSystem.thermo2();

//surfaceScalarField alphaRhoPhi2();

autoPtr<TwoFluidFoam::dragModel> pDrag = TwoFluidFoam::dragModel::New(fluidSystem);
TwoFluidFoam::dragModel& drag = pDrag();


autoPtr<saturationCurve> pSaturation = saturationCurve::New(fluidSystem);
saturationCurve& saturation = pSaturation();


/*condensationModel pCondensation = condensationModel::New
(
    fluidSystem.alpha(),
    fluidSystem.thermo2(),
    fluidSystem.U2(),
    alphaRhoPhi2,
    fluidSystem.thermo1(),
    fluidSystem.thermo2()
);*/


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

