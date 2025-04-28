{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& pMeshPhi = meshPhi.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            label ci = pFaceCells[facei];
            /*scalar ac = 0.5 * lusgsOmega * ( mag((U[ci] & pSf[facei]) - pMeshPhi[facei])  
                + a[ci]*mag(pSf[facei]));*/
            scalar ac = 0.5*lusgsOmega*(
                max(mag(U1[ci] & pSf[facei]), mag(U2[ci] & pSf[facei]))
                + max(fluid.a1()[ci]*mag(pSf[facei]), fluid.a2()[ci]*mag(pSf[facei])));

            D[ci] += ac;
        }
    }

    forAll(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rAlphaRho1  = mesh.V()[cellI]/dt.value() * dAlphaRho1[cellI];
        scalar rAlphaRho2  = mesh.V()[cellI]/dt.value() * dAlphaRho2[cellI];
        vector rAlphaRhoU1 = mesh.V()[cellI]/dt.value() * dAlphaRhoU1[cellI];
        vector rAlphaRhoU2 = mesh.V()[cellI]/dt.value() * dAlphaRhoU2[cellI];
        scalar rEpsilon1   = mesh.V()[cellI]/dt.value() * dEpsilon1[cellI];
        scalar rEpsilon2   = mesh.V()[cellI]/dt.value() * dEpsilon2[cellI];

        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index

            if (mesh.isInternalFace(faceI))
            {
                // Calculate D
                label own = owner[faceI];
                label nei = neighbour[faceI];
                //scalar dvol = mag((mesh.C()[own] - mesh.C()[nei]) & Sf[faceI]);

                scalar ac = 0.5*lusgsOmega*( 
                    max(mag(U1[cellI] & Sf[faceI]), mag(U2[cellI] & Sf[faceI])) 
                    + max(fluid.a1()[cellI] * magSf[faceI], fluid.a2()[cellI] * magSf[faceI]));

                //scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];
                scalar av = 0.0;

                D[cellI] += ac + av;

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                     // Calculate L based on A+
                    scalar ac = 0.5*lusgsOmega*( 
                        max(mag(U1[own] & Sf[faceI]), mag(U2[own] & Sf[faceI])) 
                        + max(fluid.a1()[own] * magSf[faceI], fluid.a2()[own] * magSf[faceI]));
                    //scalar av = sqr(magSf[faceI]) / dvol * nuMax[own];
                    scalar av = 0.0;
                    
                    scalar alphaRho1_1  = alphaRho1[own]  + dAlphaRho1Star[own];
                    scalar alphaRho2_1  = alphaRho2[own]  + dAlphaRho2Star[own];
                    vector alphaRhoU1_1 = alphaRhoU1[own] + dAlphaRhoU1Star[own];
                    vector alphaRhoU2_1 = alphaRhoU2[own] + dAlphaRhoU2Star[own];
                    scalar epsilon1_1   = epsilon1[own]   + dEpsilon1Star[own];
                    scalar epsilon2_1   = epsilon2[own]   + dEpsilon2Star[own];
                    
                    scalar p_1;
                    scalar alpha_1;
                    vector U1_1;
                    vector U2_1;
                    scalar T1_1;
                    scalar T2_1;

                    twoFluid::primitiveFromConservative
                    (
                        p_1,
                        alpha_1,
                        U1_1,
                        U2_1,
                        T1_1,
                        T2_1,
                        alphaRho1_1,
                        alphaRho2_1,
                        alphaRhoU1_1,
                        alphaRhoU2_1,
                        epsilon1_1,
                        epsilon2_1,
                        fluid.pInt()[celli]
                    );

                    scalar phi1_0 = U1[own] & Sf[faceI];
                    scalar phi2_0 = U2[own] & Sf[faceI];
                    scalar phi1_1 = U1_1 & Sf[faceI];
                    scalar phi2_1 = U2_1 & Sf[faceI];

                    rAlphaRho1 += (ac + av)*dAlphaRho1Star[own]
                        + 0.5*(alphaRho1_1*phi1_1 - alphaRho1[own]*phi1_0);
                    rAlphaRho2 += (ac + av)*dAlphaRho2Star[own] + 
                        0.5*(alphaRho2_1*phi2_1 - alphaRho2[own]*phi2_0);

                    rAlphaRhoU1 += (ac + av)*dAlphaRhoU1Star[own]
                        + 0.5*((alphaRhoU1_1*phi1_1 - alphaRhoU1[own]*phi1_0) + (p_1 - p[own])*Sf[faceI]);

                    rAlphaRhoU2 += (ac + av)*dAlphaRhoU2Star[own]
                        + 0.5*((alphaRhoU2_1*phi2_1 - alphaRhoU2[own]*phi2_0) + (p_1 - p[own])*Sf[faceI]);

                    rEpsilon1 += (ac + av)*dEpsilon1Star[own]
                        + 0.5*(epsilon1_1*phi1_1 + p_1*phi1_1 - (epsilon1[own]*phi1_0 + p[own]*phi1_0));

                    rEpsilon2 += (ac + av)*dEpsilon2Star[own]
                        + 0.5*(epsilon2_1*phi2_1 + p_1*phi2_1 - (epsilon1[own]*phi2_0 + p[own]*phi2_0));
                }
            }
        }

        dAlphaRho1Star[cellI]  = rAlphaRho1  / D[cellI];
        dAlphaRho2Star[cellI]  = rAlphaRho2  / D[cellI];
        dAlphaRhoU1Star[cellI] = rAlphaRhoU1 / D[cellI];
        dAlphaRhoU2Star[cellI] = rAlphaRhoU2 / D[cellI];
        dEpsilon1Star[cellI]   = rEpsilon1   / D[cellI];
        dEpsilon2Star[cellI]   = rEpsilon2   / D[cellI];
    }
}
