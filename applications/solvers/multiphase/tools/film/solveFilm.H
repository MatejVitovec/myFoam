
label wallPatchID = mesh.boundaryMesh().findPatchID("wall");

if (wallPatchID < 0)
{
    FatalErrorInFunction
        << "Patch wall not found"
        << abort(FatalError);
}


dimensionedScalar rhoLiquid
(
    "rhoLiquid",
    dimDensity,
    1000
);

const scalarField& pAlphaRhoPhi2 = twoFluidFlux.alphaRhoFlux2_pos().boundaryField()[wallPatchID];

const vectorField& Up = U1.boundaryField()[wallPatchID].patchInternalField();
const fvPatchScalarField& rhop = rho1.boundaryField()[wallPatchID];

const scalar Cf = 0.01;

forAll(h, faceI)
{
    scalar Umag = mag(Up[faceI]);

    tauG[faceI] = 0.5*rhop[faceI]*sqr(Umag)*Cf;

    Sfilm[faceI] = pAlphaRhoPhi2[faceI]/rhoLiquid.value();
}

//
dimensionedScalar muLiquid
(
    "muLiquid",
    dimViscosity,
    1e-3
);

forAll(Uf, faceI)
{
    vector n = aMesh.faceAreaNormals()[faceI];

    n /= mag(n);

    vector Ut = Up[faceI] - (Up[faceI] & n)*n;

    scalar magUt = mag(Ut);

    if (magUt > SMALL)
    {
        Uf[faceI] = ((Ut/magUt)*tauG[faceI]*h[faceI])/(2.0*muLiquid.value());
    }
    else
    {
        Uf[faceI] = Zero;
    }
}
//

edgeScalarField phiF
(
    IOobject
    (
        "phiF",
        runTime.timeName(),
        aMesh
    ),
    fac::interpolate(Uf) & aMesh.Le()
);

faScalarMatrix hEqn
(
    fam::ddt(h)
  + fam::div(phiF, h)
    ==
    Sfilm
);

hEqn.relax();

hEqn.solve();

h.max(0.0);


volScalarField hWall
(
    IOobject
    (
        "hWall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, Zero)
);

volVectorField uWall
(
    IOobject
    (
        "uWall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimLength/dimTime, vector::zero)
);

hWall.boundaryFieldRef()[wallPatchID] == h;
uWall.boundaryFieldRef()[wallPatchID] == Uf;

Info << "Film solved" << endl;