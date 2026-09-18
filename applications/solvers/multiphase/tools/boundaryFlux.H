
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
                //scalar alphaRhoPhi1Globalold = gSum(pAlphaRhoPhi1old);
                //scalar alphaRhoPhi1GlobalresSum = mag((alphaRhoPhi1Global - alphaRhoPhi1Globalold)/dt.value());

                scalar alphaRhoPhi1GlobalresSum = gSum(mag(pAlphaRhoPhi1res));
                scalar alphaRhoPhi2GlobalresSum = gSum(mag(pAlphaRhoPhi2res));

                scalar pAlphaRhoPhi1resSqr = gSum(sqr(pAlphaRhoPhi1res));
                scalar pAlphaRhoPhi2resSqr = gSum(sqr(pAlphaRhoPhi2res));
                scalar alphaRhoPhi1GlobalresNorm = Foam::sqrt(pAlphaRhoPhi1resSqr);
                scalar alphaRhoPhi2GlobalresNorm = Foam::sqrt(pAlphaRhoPhi2resSqr);
                scalar alphaRhoPhi1GlobalresNormNormalized = Foam::sqrt(pAlphaRhoPhi1resSqr/nFaces);
                scalar alphaRhoPhi2GlobalresNormNormalized = Foam::sqrt(pAlphaRhoPhi2resSqr/nFaces);

                Info << "outlet mass flux1: " << alphaRhoPhi1Global  << " , oultlet mass flux1 res: " << alphaRhoPhi1GlobalresSum << " , L2 norm1: " << alphaRhoPhi1GlobalresNorm << " , L2 norm1 normalized: " << alphaRhoPhi1GlobalresNormNormalized << endl;
                Info << "outlet mass flux2: " << alphaRhoPhi2Global  << " , oultlet mass flux2 res: " << alphaRhoPhi2GlobalresSum << " , L2 norm2: " << alphaRhoPhi2GlobalresNorm << " , L2 norm2 normalized: " << alphaRhoPhi2GlobalresNormNormalized << endl;
                if (Pstream::master() && (runTime.timeIndex() % 100) == 0)
                {
                    OFstream file
                    (
                        runTime.globalPath()/"postProcessing/boundaryOutflow.dat",
                        IOstreamOption::ASCII,
                        IOstreamOption::UNCOMPRESSED,
                        IOstreamOption::APPEND
                    );

                    file << runTime.timeIndex() << '\t' << runTime.value() << '\t' << alphaRhoPhi1Global << '\t' << alphaRhoPhi1GlobalresSum << '\t' << alphaRhoPhi1GlobalresNorm << '\t' << alphaRhoPhi1GlobalresNormNormalized
                                                                           << '\t' << alphaRhoPhi2Global << '\t' << alphaRhoPhi2GlobalresSum << '\t' << alphaRhoPhi2GlobalresNorm << '\t' << alphaRhoPhi2GlobalresNormNormalized << endl;
                }
            }
        }
    }
}
