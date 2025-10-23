//#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

TwoFluidFoam::twoFluid fluid(mesh);

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
const volScalarField& alpha2 = fluid.alpha2();
volScalarField& alpha = fluid.alpha();;
volVectorField& U1 = fluid.U1();
volVectorField& U2 = fluid.U2();
volScalarField& T1 = fluid.T1();
volScalarField& T2 = fluid.T2();

const volScalarField& pInt = fluid.pInt();

rhoThermo& thermo1 = fluid.thermo1();
rhoThermo& thermo2 = fluid.thermo2();

const volScalarField& a1 = fluid.a1();
const volScalarField& a2 = fluid.a2();

const volScalarField& e1 = thermo1.he();
const volScalarField& e2 = thermo2.he();

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

volScalarField dalpha
(
    IOobject
    (
        "dalpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dalpha", dimless, 0.0)
);

volVectorField dU1
(
    IOobject
    (
        "dU1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("dU1", dimVelocity, vector(0, 0, 0))
);

volVectorField dU2
(
    IOobject
    (
        "dU2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("dU2", dimVelocity, vector(0, 0, 0))
);


volScalarField dT1
(
    IOobject
    (
        "dT1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dT1", dimTemperature, 0.0)
);

volScalarField dT2
(
    IOobject
    (
        "dT2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("dT2", dimTemperature, 0.0)
);
