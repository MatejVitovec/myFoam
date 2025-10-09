{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), celli)
    {
        const labelList& cellFaces = mesh.cells()[celli];

        scalar rAlphaRho1  = 0;
        scalar rAlphaRho2  = 0;
        vector rAlphaRhoU1 = vector(0, 0, 0);
        vector rAlphaRhoU2 = vector(0, 0, 0);
        scalar rEpsilon1   = 0;
        scalar rEpsilon2   = 0;
            
        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                label nei = neighbour[faceI];
                label own = owner[faceI];

                if (nei != celli)
                {
                    /*const scalar ac = lusgsOmega*( 
                        max(max(mag(U1[own] & Sf[faceI]) + a1[own]*magSf[faceI],
                                mag(U2[own] & Sf[faceI]) + a2[own]*magSf[faceI]),
                            max(mag(U1[nei] & Sf[faceI]) + a1[nei]*magSf[faceI],
                                mag(U2[nei] & Sf[faceI]) + a2[nei]*magSf[faceI])));*/

                    const scalar ac = lusgsOmega*( 
                        max(mag(U1[nei] & Sf[faceI]), mag(U2[nei] & Sf[faceI])) 
                      + max(a1[nei] * magSf[faceI], a2[nei] * magSf[faceI]));

                    const scalar rhoE1 = rho1[nei]*(e1[nei] + 0.5*magSqr(U1[nei]));
                    const scalar rhoE2 = rho2[nei]*(e2[nei] + 0.5*magSqr(U2[nei]));

                    const scalar p_1     = p[nei]     + dp[nei];
                    const scalar alpha_1 = alpha[nei] + dalpha[nei];
                    const vector U1_1    = U1[nei]    + dU1[nei];
                    const vector U2_1    = U2[nei]    + dU2[nei];
                    const scalar T1_1    = T1[nei]    + dT1[nei];
                    const scalar T2_1    = T2[nei]    + dT2[nei];

                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);
                    const scalar rhoE1_1 = rho1_1*(gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1));
                    const scalar rhoE2_1 = rho2_1*(gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1));

                    const scalar phi1 = U1[nei] & Sf[faceI];
                    const scalar phi2 = U2[nei] & Sf[faceI];
                    const scalar phi1_1 = U1_1 & Sf[faceI];
                    const scalar phi2_1 = U2_1 & Sf[faceI];

                    rAlphaRho1  += 0.5*ac*(rho1_1      - rho1[nei])         + 0.5*(rho1_1*phi1_1 - rho1[nei]*phi1);
                    rAlphaRho2  += 0.5*ac*(rho2_1      - rho2[nei])         + 0.5*(rho2_1*phi2_1 - rho2[nei]*phi2); 
                    rAlphaRhoU1 += 0.5*ac*(rho1_1*U1_1 - rho1[nei]*U1[nei]) + 0.5*((rho1_1*U1_1*phi1_1 - rho1[nei]*U1[nei]*phi1) + (p_1 - p[nei])*Sf[faceI]);
                    rAlphaRhoU2 += 0.5*ac*(rho2_1*U2_1 - rho2[nei]*U2[nei]) + 0.5*((rho2_1*U2_1*phi2_1 - rho2[nei]*U2[nei]*phi2) + (p_1 - p[nei])*Sf[faceI]);
                    rEpsilon1   += 0.5*ac*(rhoE1_1     - rhoE1)             + 0.5*(rhoE1_1*phi1_1 + p_1*phi1_1 - (rhoE1*phi1 + p[nei]*phi1));
                    rEpsilon2   += 0.5*ac*(rhoE2_1     - rhoE2)             + 0.5*(rhoE2_1*phi2_1 + p_1*phi2_1 - (rhoE2*phi2 + p[nei]*phi2));
                }
            }
        }

        dp[celli]     = dpStar[celli]     + Dinv_p_alphaRho1[celli]    *rAlphaRho1 + (Dinv_p_alphaRhoU1[celli]     & rAlphaRhoU1) + Dinv_p_alphaRhoE1[celli]    *rEpsilon1 + Dinv_p_alphaRho2[celli]    *rAlphaRho2 + (Dinv_p_alphaRhoU2[celli]     & rAlphaRhoU2) + Dinv_p_alphaRhoE2[celli]    *rEpsilon2;
        dalpha[celli] = dalphaStar[celli] + Dinv_alpha_alphaRho1[celli]*rAlphaRho1 + (Dinv_alpha_alphaRhoU1[celli] & rAlphaRhoU1) + Dinv_alpha_alphaRhoE1[celli]*rEpsilon1 + Dinv_alpha_alphaRho2[celli]*rAlphaRho2 + (Dinv_alpha_alphaRhoU2[celli] & rAlphaRhoU2) + Dinv_alpha_alphaRhoE2[celli]*rEpsilon2;
        dU1[celli]    = dU1Star[celli]    + Dinv_U1_alphaRho1[celli]   *rAlphaRho1 + (Dinv_U1_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_U1_alphaRhoE1[celli]   *rEpsilon1 + Dinv_U1_alphaRho2[celli]   *rAlphaRho2 + (Dinv_U1_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_U1_alphaRhoE2[celli]   *rEpsilon2;
        dU2[celli]    = dU2Star[celli]    + Dinv_U2_alphaRho1[celli]   *rAlphaRho1 + (Dinv_U2_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_U2_alphaRhoE1[celli]   *rEpsilon1 + Dinv_U2_alphaRho2[celli]   *rAlphaRho2 + (Dinv_U2_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_U2_alphaRhoE2[celli]   *rEpsilon2;
        dT1[celli]    = dT1Star[celli]    + Dinv_T1_alphaRho1[celli]   *rAlphaRho1 + (Dinv_T1_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_T1_alphaRhoE1[celli]   *rEpsilon1 + Dinv_T1_alphaRho2[celli]   *rAlphaRho2 + (Dinv_T1_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_T1_alphaRhoE2[celli]   *rEpsilon2;
        dT2[celli]    = dT2Star[celli]    + Dinv_T2_alphaRho1[celli]   *rAlphaRho1 + (Dinv_T2_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_T2_alphaRhoE1[celli]   *rEpsilon1 + Dinv_T2_alphaRho2[celli]   *rAlphaRho2 + (Dinv_T2_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_T2_alphaRhoE2[celli]   *rEpsilon2;
    }
}
