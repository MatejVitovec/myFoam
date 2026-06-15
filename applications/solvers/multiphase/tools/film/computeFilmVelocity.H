
dimensionedScalar muLiquid
(
    "muLiquid",
    dimViscosity,
    1e-3
);

forAll(Uf, faceI)
{
    vector n =
        aMesh.faceAreaNormals()[faceI];

    n /= mag(n);

    vector Ut =
        Up[faceI]
      - (Up[faceI] & n)*n;

    scalar magUt = mag(Ut);

    if (magUt > SMALL)
    {
        Uf[faceI] =
            Ut/magUt
           *tauG[faceI]
           *h[faceI]
           /(2.0*muLiquid.value());
    }
    else
    {
        Uf[faceI] = Zero;
    }
}
