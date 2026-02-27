#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

TwoFluidFoam::twoFluid fluid(mesh);

twoFluidConservative& conservative = fluid.conservative();

twoFluidConvectiveFlux twoFluidFlux
(
    fluid.p(),
    fluid.alpha(),
    fluid.U1(),
    fluid.U2(),
    fluid.T1(),
    fluid.T2(),
    fluid.thermo1(),
    fluid.thermo2()
);


volScalarField& p = fluid.p();
const volScalarField& alpha1 = fluid.alpha1();
volScalarField& alpha2 = fluid.alpha2();
volVectorField& U1 = fluid.U1();
volVectorField& U2 = fluid.U2();
volScalarField& T1 = fluid.T1();
volScalarField& T2 = fluid.T2();

rhoThermo& thermo1 = fluid.thermo1();
rhoThermo& thermo2 = fluid.thermo2();

surfaceScalarField amaxSf
(
    "amaxSf",
    fluid.amaxSf()
);

autoPtr<TwoFluidFoam::dragModel> pDrag = TwoFluidFoam::dragModel::New(fluid);
TwoFluidFoam::dragModel& drag = pDrag();

volScalarField d
(
    IOobject
    (
        "dropletDiameter",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dropletDiameter", dimLength, 1e-5)
);


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

