
{

    forAll(mesh.boundary(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if(!isA<processorPolyPatch>(pp))
        {
            const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];

            const scalarField& pAlphaRhoPhi1 = twoFluidFlux.alphaRhoFlux1_pos().boundaryField()[patchi];
            const scalarField& pAlphaRhoPhi2 = twoFluidFlux.alphaRhoFlux2_pos().boundaryField()[patchi];

            scalar alphaRhoPhi1posGlobal = gSum(max(pAlphaRhoPhi1, scalar(0)));
            scalar alphaRhoPhi1negGlobal = gSum(min(pAlphaRhoPhi1, scalar(0)));
            scalar alphaRhoPhi2posGlobal = gSum(max(pAlphaRhoPhi2, scalar(0)));
            scalar alphaRhoPhi2negGlobal = gSum(min(pAlphaRhoPhi2, scalar(0)));
        
            Info << "patch: " << mesh.boundary()[patchi].name() << " , pos mDot.1: " << alphaRhoPhi1posGlobal << " , neg mDot.1: " << alphaRhoPhi1negGlobal << " , pos mDot.2: " << alphaRhoPhi2posGlobal << " , neg mDot.2: " << alphaRhoPhi2negGlobal << endl;
        }
    }
}
