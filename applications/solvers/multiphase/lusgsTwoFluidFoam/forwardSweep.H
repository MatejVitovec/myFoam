{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.cells(), celli)
    {
        const labelList& cellFaces = mesh.cells()[celli];

        scalar rAlphaRho1  = mesh.V()[celli]/dt.value() * rezAlphaRho1[celli];
        scalar rAlphaRho2  = mesh.V()[celli]/dt.value() * rezAlphaRho2[celli];
        vector rAlphaRhoU1 = mesh.V()[celli]/dt.value() * rezAlphaRhoU1[celli];
        vector rAlphaRhoU2 = mesh.V()[celli]/dt.value() * rezAlphaRhoU2[celli];
        scalar rEpsilon1   = mesh.V()[celli]/dt.value() * rezEpsilon1[celli];
        scalar rEpsilon2   = mesh.V()[celli]/dt.value() * rezEpsilon2[celli];

        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                if (own != celli)
                {
                    /*const scalar ac = lusgsOmega*( 
                        max(max(mag(U1[own] & Sf[faceI]) + a1[own]*magSf[faceI],
                                mag(U2[own] & Sf[faceI]) + a2[own]*magSf[faceI]),
                            max(mag(U1[nei] & Sf[faceI]) + a1[nei]*magSf[faceI],
                                mag(U2[nei] & Sf[faceI]) + a2[nei]*magSf[faceI])));*/

                    const scalar ac = lusgsOmega*(
                        max(mag(U1[own] & Sf[faceI]), mag(U2[own] & Sf[faceI])) 
                      + max(a1[own] * magSf[faceI], a2[own] * magSf[faceI]));

                    const scalar rhoE1 = rho1[own]*(e1[own] + 0.5*magSqr(U1[own]));
                    const scalar rhoE2 = rho2[own]*(e2[own] + 0.5*magSqr(U2[own]));
                    
                    const scalar p_1     = p[own]     + dpStar[own];
                    const scalar alpha_1 = alpha[own] + dalphaStar[own];
                    const vector U1_1    = U1[own]    + dU1Star[own];
                    const vector U2_1    = U2[own]    + dU2Star[own];
                    const scalar T1_1    = T1[own]    + dT1Star[own];
                    const scalar T2_1    = T2[own]    + dT2Star[own];

                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);
                    const scalar rhoE1_1 = rho1_1*(gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1));
                    const scalar rhoE2_1 = rho2_1*(gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1));

                    const scalar phi1 = U1[own] & Sf[faceI];
                    const scalar phi2 = U2[own] & Sf[faceI];
                    const scalar phi1_1 = U1_1 & Sf[faceI];
                    const scalar phi2_1 = U2_1 & Sf[faceI];

                    rAlphaRho1  += 0.5*ac*(rho1_1      - rho1[own])         + 0.5*(rho1_1*phi1_1 - rho1[own]*phi1);
                    rAlphaRho2  += 0.5*ac*(rho2_1      - rho2[own])         + 0.5*(rho2_1*phi2_1 - rho2[own]*phi2); 
                    rAlphaRhoU1 += 0.5*ac*(rho1_1*U1_1 - rho1[own]*U1[own]) + 0.5*((rho1_1*U1_1*phi1_1 - rho1[own]*U1[own]*phi1) + (p_1 - p[own])*Sf[faceI]);
                    rAlphaRhoU2 += 0.5*ac*(rho2_1*U2_1 - rho2[own]*U2[own]) + 0.5*((rho2_1*U2_1*phi2_1 - rho2[own]*U2[own]*phi2) + (p_1 - p[own])*Sf[faceI]);
                    rEpsilon1   += 0.5*ac*(rhoE1_1     - rhoE1)             + 0.5*(rhoE1_1*phi1_1 + p_1*phi1_1 - (rhoE1*phi1 + p[own]*phi1));
                    rEpsilon2   += 0.5*ac*(rhoE2_1     - rhoE2)             + 0.5*(rhoE2_1*phi2_1 + p_1*phi2_1 - (rhoE2*phi2 + p[own]*phi2));
                }
            }
        }

        dpStar[celli]     = Dinv_p_alphaRho1[celli]    *rAlphaRho1 + (Dinv_p_alphaRhoU1[celli]     & rAlphaRhoU1) + Dinv_p_alphaRhoE1[celli]    *rEpsilon1 + Dinv_p_alphaRho2[celli]    *rAlphaRho2 + (Dinv_p_alphaRhoU2[celli]     & rAlphaRhoU2) + Dinv_p_alphaRhoE2[celli]    *rEpsilon2;
        dalphaStar[celli] = Dinv_alpha_alphaRho1[celli]*rAlphaRho1 + (Dinv_alpha_alphaRhoU1[celli] & rAlphaRhoU1) + Dinv_alpha_alphaRhoE1[celli]*rEpsilon1 + Dinv_alpha_alphaRho2[celli]*rAlphaRho2 + (Dinv_alpha_alphaRhoU2[celli] & rAlphaRhoU2) + Dinv_alpha_alphaRhoE2[celli]*rEpsilon2;
        dU1Star[celli]    = Dinv_U1_alphaRho1[celli]   *rAlphaRho1 + (Dinv_U1_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_U1_alphaRhoE1[celli]   *rEpsilon1 + Dinv_U1_alphaRho2[celli]   *rAlphaRho2 + (Dinv_U1_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_U1_alphaRhoE2[celli]   *rEpsilon2;
        dU2Star[celli]    = Dinv_U2_alphaRho1[celli]   *rAlphaRho1 + (Dinv_U2_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_U2_alphaRhoE1[celli]   *rEpsilon1 + Dinv_U2_alphaRho2[celli]   *rAlphaRho2 + (Dinv_U2_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_U2_alphaRhoE2[celli]   *rEpsilon2;
        dT1Star[celli]    = Dinv_T1_alphaRho1[celli]   *rAlphaRho1 + (Dinv_T1_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_T1_alphaRhoE1[celli]   *rEpsilon1 + Dinv_T1_alphaRho2[celli]   *rAlphaRho2 + (Dinv_T1_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_T1_alphaRhoE2[celli]   *rEpsilon2;
        dT2Star[celli]    = Dinv_T2_alphaRho1[celli]   *rAlphaRho1 + (Dinv_T2_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_T2_alphaRhoE1[celli]   *rEpsilon1 + Dinv_T2_alphaRho2[celli]   *rAlphaRho2 + (Dinv_T2_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_T2_alphaRhoE2[celli]   *rEpsilon2;
    }
}
