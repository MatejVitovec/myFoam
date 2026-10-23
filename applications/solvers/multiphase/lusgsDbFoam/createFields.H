//#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

CompressibleTwoFluidFoam::compressibleTwoFluid fluid(mesh);

CompressibleTwoFluidFoam::compressibleTwoFluidConservative& conservative = fluid.conservative();

CompressibleTwoFluidFoam::compressibleTwoFluidConvectiveFlux twoFluidFlux
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

fluid.correctBoundaryCondition();
fluid.correctThermo();
fluid.correctInterfacialPressure();
fluid.correctConservative();

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

autoPtr<CompressibleTwoFluidFoam::dragModel> pDrag =
    CompressibleTwoFluidFoam::dragModel::New(fluid, d);
CompressibleTwoFluidFoam::dragModel& drag = pDrag();


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

volVectorField DU1
(
    IOobject
    (
        "DU1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::ddt(U1) + fvc::div(fvc::flux(U1), U1)
);

volVectorField DU2
(
    IOobject
    (
        "DU2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::ddt(U2) + fvc::div(fvc::flux(U2), U2)
);
