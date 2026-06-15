
{
    massSourceWall2 = dimensionedScalar("massSourceWall2zero", dimMass/dimTime/dimVolume, 0.0);

    forAll(mesh.boundary(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if(isA<wallPolyPatch>(pp))
        {
            const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];

            forAll(pFaceCells, i)
            {
                label cellI = pFaceCells[i];

                const scalar volAlphaRhoPhi2 = alpha2[cellI]*rho2[cellI]*(U2[cellI] & pSf[i]);

                if (volAlphaRhoPhi2 > 0)
                {
                    massSourceWall2[cellI] += volAlphaRhoPhi2/mesh.V()[cellI];
                }
            }

            Info << "phase2 outflow sink: " << gSum(massSourceWall2) << endl;
        }
    }
}