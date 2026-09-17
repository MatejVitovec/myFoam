
{

    forAll(mesh.boundary(), patchi)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        if(!isA<processorPolyPatch>(pp))
        {
            const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
            const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
            const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];

            //const scalarField& pAlphaRhoPhi1 = twoFluidFlux.alphaRhoFlux1_pos().boundaryField()[patchi];
            const scalarField& pAlphaRhoPhi2 = twoFluidFlux.alphaRhoFlux2_pos().boundaryField()[patchi];

            //scalar alphaRhoPhi1posGlobal = gSum(max(pAlphaRhoPhi1, scalar(0)));
            //scalar alphaRhoPhi1negGlobal = gSum(min(pAlphaRhoPhi1, scalar(0)));
            //scalar alphaRhoPhi2posGlobal = gSum(max(pAlphaRhoPhi2, scalar(0)));
            //scalar alphaRhoPhi2negGlobal = gSum(min(pAlphaRhoPhi2, scalar(0)));

            //Info << "patch: " << mesh.boundary()[patchi].name() << " , pos mDot.1: " << alphaRhoPhi1posGlobal << " , neg mDot.1: " << alphaRhoPhi1negGlobal << " , pos mDot.2: " << alphaRhoPhi2posGlobal << " , neg mDot.2: " << alphaRhoPhi2negGlobal << endl;

            if(pp.name() == "outlet")
            {
                const scalarField& pAlphaRhoPhi1    = twoFluidFlux.alphaRhoFlux1_pos().boundaryField()[patchi];
                const scalarField& pAlphaRhoPhi1old = twoFluidFlux.alphaRhoFlux1_pos().oldTime().boundaryField()[patchi];
                const scalarField& pAlphaRhoPhi2    = twoFluidFlux.alphaRhoFlux2_pos().boundaryField()[patchi];
                const scalarField& pAlphaRhoPhi2old = twoFluidFlux.alphaRhoFlux2_pos().oldTime().boundaryField()[patchi];

                const scalarField pAlphaRhoPhi1res = (pAlphaRhoPhi1 - pAlphaRhoPhi1old)/dt.value();
                const scalarField pAlphaRhoPhi2res = (pAlphaRhoPhi2 - pAlphaRhoPhi2old)/dt.value();
                const label nFaces = returnReduce(mesh.boundary()[pp.name()].size(), sumOp<label>());

                scalar alphaRhoPhi1Global = gSum(pAlphaRhoPhi1);
                scalar alphaRhoPhi2Global = gSum(pAlphaRhoPhi2);
                scalar alphaRhoPhi1Globalold = gSum(pAlphaRhoPhi1old);
                //scalar alphaRhoPhi1GlobalresSum = mag((alphaRhoPhi1Global - alphaRhoPhi1Globalold)/dt.value());

                scalar alphaRhoPhi1SfGlobalresSum = gSum(mag(pAlphaRhoPhi1res*pMagSf));
                scalar alphaRhoPhi2SfGlobalresSum = gSum(mag(pAlphaRhoPhi2res*pMagSf));

                scalar pAlphaRhoPhi1resSfSqr = gSum(sqr(pAlphaRhoPhi1res*pMagSf));
                scalar pAlphaRhoPhi2resSfSqr = gSum(sqr(pAlphaRhoPhi2res*pMagSf));
                scalar alphaRhoPhi1SfGlobalresNorm = Foam::sqrt(pAlphaRhoPhi1resSfSqr);
                scalar alphaRhoPhi2SfGlobalresNorm = Foam::sqrt(pAlphaRhoPhi2resSfSqr);
                scalar alphaRhoPhi1SfGlobalresNormNormalized = Foam::sqrt(pAlphaRhoPhi1resSfSqr/nFaces);
                scalar alphaRhoPhi2SfGlobalresNormNormalized = Foam::sqrt(pAlphaRhoPhi2resSfSqr/nFaces);

                Info << "outlet mass flux1: " << alphaRhoPhi1Global  << " , oultlet mass flux1 res: " << alphaRhoPhi1SfGlobalresSum << " , L2 norm1: " << alphaRhoPhi1SfGlobalresNorm << " , L2 norm1 normalized: " << alphaRhoPhi1SfGlobalresNormNormalized << endl;
                Info << "outlet mass flux2: " << alphaRhoPhi2Global  << " , oultlet mass flux2 res: " << alphaRhoPhi2SfGlobalresSum << " , L2 norm2: " << alphaRhoPhi2SfGlobalresNorm << " , L2 norm2 normalized: " << alphaRhoPhi2SfGlobalresNormNormalized << endl;
                if (Pstream::master() && (runTime.timeIndex() % 100) == 0)
                {
                    OFstream file
                    (
                        runTime.globalPath()/"postProcessing/boundaryOutflow.dat",
                        IOstreamOption::ASCII,
                        IOstreamOption::UNCOMPRESSED,
                        IOstreamOption::APPEND
                    );

                    file << runTime.timeIndex() << '\t' << runTime.value() << '\t' << alphaRhoPhi1Global << '\t' << alphaRhoPhi1SfGlobalresSum << '\t' << alphaRhoPhi1SfGlobalresNorm << '\t' << alphaRhoPhi1SfGlobalresNormNormalized
                                                                           << '\t' << alphaRhoPhi2Global << '\t' << alphaRhoPhi2SfGlobalresSum << '\t' << alphaRhoPhi2SfGlobalresNorm << '\t' << alphaRhoPhi2SfGlobalresNormNormalized << endl;
                }
            }
        }
    }
}
