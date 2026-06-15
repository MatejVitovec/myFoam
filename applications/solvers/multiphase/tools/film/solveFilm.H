
label wallPatchID = mesh.boundaryMesh().findPatchID("wall");

if (wallPatchID < 0)
{
    FatalErrorInFunction
        << "Patch wall not found"
        << abort(FatalError);
}

const fvPatchVectorField& Up = U1.boundaryField()[wallPatchID];
const fvPatchScalarField& rhop = rho1.boundaryField()[wallPatchID];

const scalar Cf = 0.005;

forAll(h, faceI)
{
    scalar Umag = mag(Up[faceI]);

    tauG[faceI] = 0.5*rhop[faceI]*sqr(Umag)*Cf;
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
