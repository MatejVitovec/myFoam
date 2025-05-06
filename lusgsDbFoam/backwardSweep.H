
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rAlphaRho1  = 0;
        scalar rAlphaRho2  = 0;
        vector rAlphaRhoU1 = vector(0, 0, 0);
        vector rAlphaRhoU2 = vector(0, 0, 0);
        scalar rEpsilon1   = 0;
        scalar rEpsilon2   = 0;
            
        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index
            if (mesh.isInternalFace(faceI))
            {
                label nei = neighbour[faceI];
                label own = owner[faceI];

                // we use the mesh ordering property: own < nei
                if (nei != cellI)  // then data[nei] belongs to upper triangle
                {
                    //scalar dvol = mag( (mesh.C()[own] - mesh.C()[nei]) & Sf[faceI] );
                    /*const scalar ac = 0.5*lusgsOmega*( 
                        max(mag(U1[nei] & Sf[faceI]), mag(U2[nei] & Sf[faceI])) 
                      + max(a1[nei] * magSf[faceI], a2[nei] * magSf[faceI]));*/
                    const scalar ac = 0.5*lusgsOmega*( 
                        max(max(mag(U1[own] & Sf[faceI]) + a1[own]*magSf[faceI],
                                mag(U2[own] & Sf[faceI]) + a2[own]*magSf[faceI]),
                            max(mag(U1[nei] & Sf[faceI]) + a1[nei]*magSf[faceI],
                                mag(U2[nei] & Sf[faceI]) + a2[nei]*magSf[faceI])));
                    //scalar av = sqr(magSf[faceI]) / dvol * nuMax[nei];
                    const scalar av = 0.0;

                    const scalar alphaRho1_1  = alphaRho1[nei]  + dAlphaRho1[nei];
                    const scalar alphaRho2_1  = alphaRho2[nei]  + dAlphaRho2[nei];
                    const vector alphaRhoU1_1 = alphaRhoU1[nei] + dAlphaRhoU1[nei];
                    const vector alphaRhoU2_1 = alphaRhoU2[nei] + dAlphaRhoU2[nei];
                    const scalar epsilon1_1   = epsilon1[nei]   + dEpsilon1[nei];
                    const scalar epsilon2_1   = epsilon2[nei]   + dEpsilon2[nei];

                    /*scalar fluxAlphaRho1_1;
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
                        fluid.pInt()[nei],
                        Sf[faceI],
                        p[nei],
                        T1[nei],
                        T2[nei]
                    );

                    const scalar alphaRhoE1_0 = epsilon1[nei] - alpha[nei]*fluid.pInt()[nei];
                    const scalar alphaRhoE2_0 = epsilon2[nei] - (1.0 - alpha[nei])*fluid.pInt()[nei];
                    const scalar alphaP1_0 = alpha[nei]*p[nei];
                    const scalar alphaP2_0 = (1.0 - alpha[nei])*p[nei];
                    const scalar phi1_0 = U1[nei] & Sf[faceI];
                    const scalar phi2_0 = U2[nei] & Sf[faceI];

                    rAlphaRho1 += (ac + av)*dAlphaRho1Star[nei]
                        + 0.5*(fluxAlphaRho1_1 - alphaRho1[nei]*phi1_0);
                    rAlphaRho2 += (ac + av)*dAlphaRho2Star[nei]
                        + 0.5*(fluxAlphaRho2_1 - alphaRho2[nei]*phi2_0);

                    rAlphaRhoU1 += (ac + av)*dAlphaRhoU1Star[nei]
                        + 0.5*(fluxAlphaRhoU1_1 - (alphaRhoU1[nei]*phi1_0 + alphaP1_0*Sf[faceI]));
                    rAlphaRhoU2 += (ac + av)*dAlphaRhoU2Star[nei]
                        + 0.5*(fluxAlphaRhoU2_1 - (alphaRhoU2[nei]*phi2_0 + alphaP2_0*Sf[faceI]));

                    rEpsilon1 += (ac + av)*dEpsilon1Star[nei]
                        + 0.5*(fluxAlphaRhoE1_1 - (alphaRhoE1_0 + alphaP1_0)*phi1_0);
                    rEpsilon2 += (ac + av)*dEpsilon2Star[nei]
                        + 0.5*(fluxAlphaRhoE2_1 - (alphaRhoE2_0 + alphaP2_0)*phi2_0);*/


                    //this also serves as an init value for Newton iteration in primitiveFromConservative function
                    scalar p_1 = p[nei];
                    scalar T1_1 = T1[nei];
                    scalar T2_1 = T2[nei];

                    scalar alpha1_1;
                    vector U1_1;
                    vector U2_1;

                    const scalar pInt = fluid.pInt()[nei];

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

                    //fluid.blendVanishingFluid(alpha1_1, U1_1, U2_1, T1_1, T2_1);

                    const scalar alphaRhoE1_0 = epsilon1[nei] - alpha[nei]*pInt;
                    const scalar alphaRhoE2_0 = epsilon2[nei] - (1.0 - alpha[nei])*pInt;
                    const scalar alphaRhoE1_1 = epsilon1_1 - alpha1_1*pInt;
                    const scalar alphaRhoE2_1 = epsilon2_1 - (1.0 - alpha1_1)*pInt;

                    const scalar alphaP1_0 = alpha[nei]*p[nei];
                    const scalar alphaP2_0 = (1.0 - alpha[nei])*p[nei];
                    const scalar alphaP1_1 = alpha1_1*p_1;
                    const scalar alphaP2_1 = (1.0 - alpha1_1)*p_1;

                    const scalar phi1_0 = U1[nei] & Sf[faceI];
                    const scalar phi2_0 = U2[nei] & Sf[faceI];
                    const scalar phi1_1 = U1_1 & Sf[faceI];
                    const scalar phi2_1 = U2_1 & Sf[faceI];

                    rAlphaRho1 += (ac + av) * dAlphaRho1[nei] 
                        - 0.5*(alphaRho1_1*phi1_1 - alphaRho1[nei]*phi1_0);
                    rAlphaRho2 += (ac + av) * dAlphaRho2[nei] 
                        - 0.5*(alphaRho2_1*phi2_1 - alphaRho2[nei]*phi2_0);

                    rAlphaRhoU1 += (ac + av) * dAlphaRhoU1[nei]
                        - 0.5*((alphaRhoU1_1*phi1_1 - alphaRhoU1[nei]*phi1_0) + (alphaP1_1 - alphaP1_0)*Sf[faceI]);
                    rAlphaRhoU2 += (ac + av) * dAlphaRhoU2[nei]
                        - 0.5*((alphaRhoU2_1*phi2_1 - alphaRhoU2[nei]*phi2_0) + (alphaP2_1 - alphaP2_0)*Sf[faceI]);

                    rEpsilon1 += (ac + av) * dEpsilon1[nei]
                        - 0.5*((alphaRhoE1_1 + alphaP1_1)*phi1_1 - (alphaRhoE1_0 + alphaP1_0)*phi1_0);
                    rEpsilon2 += (ac + av) * dEpsilon2[nei]
                        - 0.5*((alphaRhoE2_1 + alphaP2_1)*phi2_1 - (alphaRhoE2_0 + alphaP2_0)*phi2_0);
                }
            }
        }
        
        dAlphaRho1[cellI]  = dAlphaRho1Star[cellI]  + rAlphaRho1  / D[cellI];
        dAlphaRho2[cellI]  = dAlphaRho2Star[cellI]  + rAlphaRho2  / D[cellI];
        dAlphaRhoU1[cellI] = dAlphaRhoU1Star[cellI] + rAlphaRhoU1 / D[cellI];
        dAlphaRhoU2[cellI] = dAlphaRhoU2Star[cellI] + rAlphaRhoU2 / D[cellI];
        dEpsilon1[cellI]   = dEpsilon1Star[cellI]   + rEpsilon1   / D[cellI];
        dEpsilon2[cellI]   = dEpsilon2Star[cellI]   + rEpsilon2   / D[cellI];
    }
}
