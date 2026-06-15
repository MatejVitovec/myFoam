
label wallPatchID =
    mesh.boundaryMesh().findPatchID("wall");

if (wallPatchID < 0)
{
    FatalErrorInFunction
        << "Patch wall not found"
        << abort(FatalError);
}

const fvPatchVectorField& Up =
    U1.boundaryField()[wallPatchID];

const fvPatchScalarField& rhop =
    rho1.boundaryField()[wallPatchID];

scalar Cf = 0.005;

forAll(h, faceI)
{
    scalar Umag = mag(Up[faceI]);

    tauG[faceI] =
        0.5
       *rhop[faceI]
       *sqr(Umag)
       *Cf;
}
