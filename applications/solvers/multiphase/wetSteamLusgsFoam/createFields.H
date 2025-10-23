//#include "createRDeltaT.H"

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

volScalarField rho
(
    "rho",
    fluid.rho()
);

surfaceScalarField rhoMPhi2
(
    "rhoMPhi2",
    linearInterpolate(rho*U2) & mesh.Sf()
);

autoPtr<TwoFluidFoam::dragModel> pDrag = TwoFluidFoam::dragModel::New(fluid);
TwoFluidFoam::dragModel& drag = pDrag();

autoPtr<Foam::WetSteam::condensationModel> pCondensation = Foam::WetSteam::condensationModel::New
(
    fluid.alpha2(),
    fluid.thermo2().rho(),
    fluid.U2(),
    rhoMPhi2,
    fluid.thermo1(),
    fluid.thermo2(),
    fluid
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

/*volVectorField DU1
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
);*/
