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

const volScalarField& rho1 = thermo1.rho();
const volScalarField& rho2 = thermo2.rho();

surfaceScalarField amaxSf
(
    "amaxSf",
    fluid.amaxSf()
);

const volScalarField alphaM
(
    IOobject
    (
        "alphaMixture",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaMixtureIsOne", dimless, 1.0)
);

volScalarField rhoM
(
    "rhoMixture",
    fluid.rho()
);

surfaceScalarField rhoMPhi2
(
    "rhoMPhi2",
    linearInterpolate(rhoM*U2) & mesh.Sf()
);

autoPtr<TwoFluidFoam::dragModel> pDrag = TwoFluidFoam::dragModel::New(fluid);
TwoFluidFoam::dragModel& drag = pDrag();

saturation satur = saturation(fluid.subDict("condensation"), thermo1);

autoPtr<Foam::WetSteam::condensationModel> pCondensation = Foam::WetSteam::condensationModel::New
(
    fluid.alpha2(),
    alphaM,
    rhoM,
    fluid.U2(),
    rhoMPhi2,
    fluid.thermo1(),
    fluid.thermo2(),
    satur,
    fluid
);
Foam::WetSteam::condensationModel& condensation = pCondensation();

/*volScalarField d
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
);*/

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
