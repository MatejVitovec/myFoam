
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rRho  = 0;
        vector rRhoU = vector(0,0,0);
        scalar rRhoE = 0;

        scalar rAlphaRho1  = 0;
        scalar rAlphaRho1  = 0;
        vector rAlphaRhoU1 = vector(0,0,0);
        vector rAlphaRhoU2 = vector(0,0,0);
        scalar rEpsilon1   = 0;
        scalar rEpsilon1   = 0;
            
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
                    scalar ac = 0.5*lusgsOmega*( 
                        max(mag(U1[nei] & Sf[faceI]), mag(U2[nei] & Sf[faceI])) 
                      + max(a1[nei] * magSf[faceI], a2[nei] * magSf[faceI]));
                    //scalar av = sqr(magSf[faceI]) / dvol * nuMax[nei];
                    scalar av = 0.0;

                    scalar alphaRho1_1  = alphaRho1[nei]  + dAlphaRho1[nei];
                    scalar alphaRho2_1  = alphaRho2[nei]  + dAlphaRho2[nei];
                    vector alphaRhoU1_1 = alphaRhoU1[nei] + dAlphaRhoU1[nei];
                    vector alphaRhoU2_1 = alphaRhoU2[nei] + dAlphaRhoU2[nei];
                    scalar epsilon1_1   = epsilon1[nei]   + dEpsilon1[nei];
                    scalar epsilon2_1   = epsilon2[nei]   + dEpsilon2[nei];

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

                    rAlphaRho1 += (ac + av) * dAlphaRho1[nei] 
                        - 0.5*(alphaRho1_1*phi1_1 - alphaRho1[nei]*ph1_0);
                    rAlphaRho2 += (ac + av) * dAlphaRho2[nei] 
                        - 0.5*(alphaRho2_1*phi2_1 - alphaRho2[nei]*ph2_0);

                    rAlphaRhoU1 += (ac + av) * dAlphaRhoU1[nei]
                        - 0.5*((alphaRhoU1_1*phi1_1 - alphaRhoU1[nei]*phi1_0) + (p_1 - p[nei])*Sf[faceI]);
                    rAlphaRhoU2 += (ac + av) * dAlphaRhoU2[nei]
                        - 0.5*((alphaRhoU2_1*phi2_1 - alphaRhoU2[nei]*phi2_0) + (p_1 - p[nei])*Sf[faceI]);

                    rEpsilon1 += (ac + av) * dEpsilon1[nei]
                        - 0.5*((epsilon1_1*phi1_1 + p_1*phi1_1) - (epsilon1[nei]*phi1_0 + p[nei]*phi1_0));
                    rEpsilon2 += (ac + av) * dEpsilon2[nei]
                        - 0.5*((epsilon2_1*phi2_1 + p_1*phi2_1) - (epsilon2[nei]*phi2_0 + p[nei]*phi2_0));
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
