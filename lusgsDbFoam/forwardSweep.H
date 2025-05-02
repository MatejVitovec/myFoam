{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            label ci = pFaceCells[facei];
            /*scalar ac = 0.5 * lusgsOmega * ( mag((U[ci] & pSf[facei]) - pMeshPhi[facei])  
                + a[ci]*mag(pSf[facei]));*/
            scalar ac = 0.5*lusgsOmega*(
                max(mag(U1[ci] & pSf[facei]), mag(U2[ci] & pSf[facei]))
                + max(a1[ci]*mag(pSf[facei]), a2[ci]*mag(pSf[facei])));

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

                /*const scalar ac = 0.5*lusgsOmega*( 
                    max(mag(U1[cellI] & Sf[faceI]), mag(U2[cellI] & Sf[faceI])) 
                    + max(a1[cellI]*magSf[faceI], a2[cellI]*magSf[faceI]));*/
                const scalar ac = 0.5*lusgsOmega*( 
                    max(max(mag(U1[cellI] & Sf[faceI]) + a1[cellI]*magSf[faceI],
                            mag(U2[cellI] & Sf[faceI]) + a2[cellI]*magSf[faceI]),
                        max(mag(U1[nei] & Sf[faceI]) + a1[nei]*magSf[faceI],
                            mag(U2[nei] & Sf[faceI]) + a2[nei]*magSf[faceI])));

                //scalar av = sqr(magSf[faceI]) / dvol * nuMax[cellI];
                const scalar av = 0.0;

                D[cellI] += ac + av;

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                     // Calculate L based on A+
                    /*const scalar ac = 0.5*lusgsOmega*( 
                        max(mag(U1[own] & Sf[faceI]), mag(U2[own] & Sf[faceI])) 
                        + max(a1[own] * magSf[faceI], a2[own] * magSf[faceI]));*/

                    //scalar av = sqr(magSf[faceI]) / dvol * nuMax[own];
                    //const scalar av = 0.0;
                    
                    const scalar alphaRho1_1  = alphaRho1[own]  + dAlphaRho1Star[own];
                    const scalar alphaRho2_1  = alphaRho2[own]  + dAlphaRho2Star[own];
                    const vector alphaRhoU1_1 = alphaRhoU1[own] + dAlphaRhoU1Star[own];
                    const vector alphaRhoU2_1 = alphaRhoU2[own] + dAlphaRhoU2Star[own];
                    const scalar epsilon1_1   = epsilon1[own]   + dEpsilon1Star[own];
                    const scalar epsilon2_1   = epsilon2[own]   + dEpsilon2Star[own];

                    scalar fluxAlphaRho1_1;
                    scalar fluxAlphaRho2_1;
                    vector fluxAlphaRhoU1_1;
                    vector fluxAlphaRhoU2_1;
                    scalar fluxAlphaRhoE1_1;
                    scalar fluxAlphaRhoE2_1;

                    fluid.fluxFromConservative
                    (
                        fluxAlphaRho1_1,
                        fluxAlphaRho2_1,
                        fluxAlphaRhoU1_1,
                        fluxAlphaRhoU2_1,
                        fluxAlphaRhoE1_1,
                        fluxAlphaRhoE2_1,
                        alphaRho1_1,
                        alphaRho2_1,
                        alphaRhoU1_1,
                        alphaRhoU2_1,
                        epsilon1_1,
                        epsilon2_1,
                        fluid.pInt()[own],
                        Sf[faceI],
                        p[own],
                        T1[own],
                        T2[own]
                    );

                    const scalar alphaRhoE1_0 = epsilon1[own] - alpha[own]*fluid.pInt()[own];
                    const scalar alphaRhoE2_0 = epsilon2[own] - (1.0 - alpha[own])*fluid.pInt()[own];
                    const scalar alphaP1_0 = alpha[own]*p[own];
                    const scalar alphaP2_0 = (1.0 - alpha[own])*p[own];
                    const scalar phi1_0 = U1[own] & Sf[faceI];
                    const scalar phi2_0 = U2[own] & Sf[faceI];

                    rAlphaRho1 += (ac + av)*dAlphaRho1Star[own]
                        + 0.5*(fluxAlphaRho1_1 - alphaRho1[own]*phi1_0);
                    rAlphaRho2 += (ac + av)*dAlphaRho2Star[own]
                        + 0.5*(fluxAlphaRho2_1 - alphaRho2[own]*phi2_0);

                    rAlphaRhoU1 += (ac + av)*dAlphaRhoU1Star[own]
                        + 0.5*(fluxAlphaRhoU1_1 - (alphaRhoU1[own]*phi1_0 + alphaP1_0*Sf[faceI]));
                    rAlphaRhoU2 += (ac + av)*dAlphaRhoU2Star[own]
                        + 0.5*(fluxAlphaRhoU2_1 - (alphaRhoU2[own]*phi2_0 + alphaP2_0*Sf[faceI]));

                    rEpsilon1 += (ac + av)*dEpsilon1Star[own]
                        + 0.5*(fluxAlphaRhoE1_1 - (alphaRhoE1_0 + alphaP1_0)*phi1_0);
                    rEpsilon2 += (ac + av)*dEpsilon2Star[own]
                        + 0.5*(fluxAlphaRhoE2_1 - (alphaRhoE2_0 + alphaP2_0)*phi2_0);

                    /*
                    //this also serves as an init value for Newton iteration in primitiveFromConservative function
                    scalar p_1 = p[own];
                    scalar T1_1 = T1[own];
                    scalar T2_1 = T2[own];

                    scalar alpha1_1;
                    vector U1_1;
                    vector U2_1;

                    const scalar pInt = fluid.pInt()[own];

                    fluid.primitiveFromConservative
                    (
                        p_1,
                        alpha1_1,
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
                        pInt
                    );

                    fluid.blendVanishingFluid(alpha1_1, U1_1, U2_1, T1_1, T2_1); //TODO i pro conservative

                    const scalar alphaRhoE1_0 = epsilon1[own] - alpha[own]*pInt;
                    const scalar alphaRhoE2_0 = epsilon2[own] - (1.0 - alpha[own])*pInt;
                    const scalar alphaRhoE1_1 = epsilon1_1 - alpha1_1*pInt;
                    const scalar alphaRhoE2_1 = epsilon2_1 - (1.0 - alpha1_1)*pInt;

                    const scalar alphaP1_0 = alpha[own]*p[own];
                    const scalar alphaP2_0 = (1.0 - alpha[own])*p[own];
                    const scalar alphaP1_1 = alpha1_1*p_1;
                    const scalar alphaP2_1 = (1.0 - alpha1_1)*p_1;

                    const scalar phi1_0 = U1[own] & Sf[faceI];
                    const scalar phi2_0 = U2[own] & Sf[faceI];
                    const scalar phi1_1 = U1_1 & Sf[faceI];
                    const scalar phi2_1 = U2_1 & Sf[faceI];

                    rAlphaRho1 += (ac + av)*dAlphaRho1Star[own]
                        + 0.5*(alphaRho1_1*phi1_1 - alphaRho1[own]*phi1_0);
                    rAlphaRho2 += (ac + av)*dAlphaRho2Star[own]
                        + 0.5*(alphaRho2_1*phi2_1 - alphaRho2[own]*phi2_0);

                    rAlphaRhoU1 += (ac + av)*dAlphaRhoU1Star[own]
                        + 0.5*((alphaRhoU1_1*phi1_1 - alphaRhoU1[own]*phi1_0) + (alphaP1_1 - alphaP1_0)*Sf[faceI]);
                    rAlphaRhoU2 += (ac + av)*dAlphaRhoU2Star[own]
                        + 0.5*((alphaRhoU2_1*phi2_1 - alphaRhoU2[own]*phi2_0) + (alphaP2_1 - alphaP2_0)*Sf[faceI]);

                    rEpsilon1 += (ac + av)*dEpsilon1Star[own]
                        + 0.5*((alphaRhoE1_1 + alphaP1_1)*phi1_1 - (alphaRhoE1_0 + alphaP1_0)*phi1_0);
                    rEpsilon2 += (ac + av)*dEpsilon2Star[own]
                        + 0.5*((alphaRhoE2_1 + alphaP2_1)*phi2_1 - (alphaRhoE2_0 + alphaP2_0)*phi2_0);*/
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
