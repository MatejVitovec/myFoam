
Info<< "Reading finite-area mesh" << endl;

label patchID = mesh.boundaryMesh().findPatchID("wall");

if (patchID < 0)
{
    FatalErrorInFunction
        << "Patch wall not found"
        << abort(FatalError);
}

const polyPatch& pp = mesh.boundaryMesh()[patchID];

faMesh aMesh(pp);

Info<< "Creating film thickness field h" << endl;

areaScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        aMesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedScalar("h0", dimLength, 1e-20)
);

Info<< "Creating film velocity field" << endl;

areaVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        aMesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh,
    dimensionedVector("zero", dimVelocity, Zero)
);
