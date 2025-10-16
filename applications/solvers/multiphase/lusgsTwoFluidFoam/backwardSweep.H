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
            label facei = cellFaces[i];

            if (mesh.isInternalFace(facei))
            {
                label nei = neighbour[facei];
                //label own = owner[facei];

                if (nei != celli)
                {
                    /*const scalar ac = lusgsOmega*( 
                        max(max(mag(U1[own] & Sf[facei]) + a1[own]*magSf[facei],
                                mag(U2[own] & Sf[facei]) + a2[own]*magSf[facei]),
                            max(mag(U1[nei] & Sf[facei]) + a1[nei]*magSf[facei],
                                mag(U2[nei] & Sf[facei]) + a2[nei]*magSf[facei])));*/

                    //const scalar ac = lusgsOmega*(max(mag(U1[nei] & Sf[facei]), mag(U2[nei] & Sf[facei])) + max(a1[nei] * magSf[facei], a2[nei] * magSf[facei]));

                    const scalar ac1 = lusgsOmega*(mag(U1[nei] & Sf[facei]) + a1[nei]*magSf[facei]);
                    const scalar ac2 = lusgsOmega*(mag(U2[nei] & Sf[facei]) + a2[nei]*magSf[facei]);

                    const scalar p_1 = p[nei]     + dp[nei];
                    scalar alpha_1   = alpha[nei] + dalpha[nei];
                    vector U1_1      = U1[nei]    + dU1[nei];
                    vector U2_1      = U2[nei]    + dU2[nei];
                    scalar T1_1      = T1[nei]    + dT1[nei];
                    scalar T2_1      = T2[nei]    + dT2[nei];

                    //fluid.blendVanishingFluid(alpha_1, U1_1, U2_1, T1_1, T2_1);

                    const scalar alpha1_1 = 1.0 - alpha_1;
                    const scalar alpha2_1 = alpha_1;
                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);

                    const scalar E1   = e1[nei]                 + 0.5*magSqr(U1[nei]);
                    const scalar E2   = e2[nei]                 + 0.5*magSqr(U2[nei]);
                    const scalar E1_1 = gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1);
                    const scalar E2_1 = gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1);

                    /*const scalar rhoE1   = rho1[nei]*(e1[nei]              + 0.5*magSqr(U1[nei]));
                    const scalar rhoE2   = rho2[nei]*(e2[nei]              + 0.5*magSqr(U2[nei]));
                    const scalar rhoE1_1 = rho1_1*(gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1));
                    const scalar rhoE2_1 = rho2_1*(gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1));*/

                    const scalar alphaRhoPhi1   = alpha1[nei]*rho1[nei]*(U1[nei] & Sf[facei]);
                    const scalar alphaRhoPhi2   = alpha2[nei]*rho2[nei]*(U2[nei] & Sf[facei]);
                    const scalar alphaRhoPhi1_1 = alpha1_1*rho1_1*(U1_1 & Sf[facei]);
                    const scalar alphaRhoPhi2_1 = alpha2_1*rho2_1*(U2_1 & Sf[facei]);

                    rAlphaRho1  += 0.5*ac1*(alpha1_1*rho1_1                    - alpha1[nei]*rho1[nei])                  + 0.5*(alphaRhoPhi1_1                                 - alphaRhoPhi1);
                    rAlphaRho2  += 0.5*ac2*(alpha2_1*rho2_1                    - alpha2[nei]*rho2[nei])                  + 0.5*(alphaRhoPhi2_1                                 - alphaRhoPhi2); 
                    rAlphaRhoU1 += 0.5*ac1*(alpha1_1*rho1_1*U1_1               - alpha1[nei]*rho1[nei]*U1[nei])          + 0.5*((alphaRhoPhi1_1*U1_1 + alpha1_1*p_1*Sf[facei]) - (alphaRhoPhi1*U1[nei] + alpha1[nei]*p[nei]*Sf[facei]));
                    rAlphaRhoU2 += 0.5*ac2*(alpha2_1*rho2_1*U2_1               - alpha2[nei]*rho2[nei]*U2[nei])          + 0.5*((alphaRhoPhi2_1*U2_1 + alpha2_1*p_1*Sf[facei]) - (alphaRhoPhi2*U2[nei] + alpha2[nei]*p[nei]*Sf[facei]));
                    rEpsilon1   += 0.5*ac1*(alpha1_1*(rho1_1*E1_1 + pInt[nei]) - alpha1[nei]*(rho1[nei]*E1 + pInt[nei])) + 0.5*(alphaRhoPhi1_1*(E1_1 + p_1/rho1_1)             - alphaRhoPhi1*(E1 + p[nei]/rho1[nei]));
                    rEpsilon2   += 0.5*ac2*(alpha2_1*(rho2_1*E2_1 + pInt[nei]) - alpha2[nei]*(rho2[nei]*E2 + pInt[nei])) + 0.5*(alphaRhoPhi2_1*(E2_1 + p_1/rho2_1)             - alphaRhoPhi2*(E2 + p[nei]/rho2[nei]));
                }
            }
        }

        /*dp[celli]     = dpStar[celli]     + Dinv_p_alphaRho1[celli]    *rAlphaRho1 + (Dinv_p_alphaRhoU1[celli]     & rAlphaRhoU1) + Dinv_p_alphaRhoE1[celli]    *rEpsilon1 + Dinv_p_alphaRho2[celli]    *rAlphaRho2 + (Dinv_p_alphaRhoU2[celli]     & rAlphaRhoU2) + Dinv_p_alphaRhoE2[celli]    *rEpsilon2;
        dalpha[celli] = dalphaStar[celli] + Dinv_alpha_alphaRho1[celli]*rAlphaRho1 + (Dinv_alpha_alphaRhoU1[celli] & rAlphaRhoU1) + Dinv_alpha_alphaRhoE1[celli]*rEpsilon1 + Dinv_alpha_alphaRho2[celli]*rAlphaRho2 + (Dinv_alpha_alphaRhoU2[celli] & rAlphaRhoU2) + Dinv_alpha_alphaRhoE2[celli]*rEpsilon2;
        dU1[celli]    = dU1Star[celli]    + Dinv_U1_alphaRho1[celli]   *rAlphaRho1 + (Dinv_U1_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_U1_alphaRhoE1[celli]   *rEpsilon1 + Dinv_U1_alphaRho2[celli]   *rAlphaRho2 + (Dinv_U1_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_U1_alphaRhoE2[celli]   *rEpsilon2;
        dU2[celli]    = dU2Star[celli]    + Dinv_U2_alphaRho1[celli]   *rAlphaRho1 + (Dinv_U2_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_U2_alphaRhoE1[celli]   *rEpsilon1 + Dinv_U2_alphaRho2[celli]   *rAlphaRho2 + (Dinv_U2_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_U2_alphaRhoE2[celli]   *rEpsilon2;
        dT1[celli]    = dT1Star[celli]    + Dinv_T1_alphaRho1[celli]   *rAlphaRho1 + (Dinv_T1_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_T1_alphaRhoE1[celli]   *rEpsilon1 + Dinv_T1_alphaRho2[celli]   *rAlphaRho2 + (Dinv_T1_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_T1_alphaRhoE2[celli]   *rEpsilon2;
        dT2[celli]    = dT2Star[celli]    + Dinv_T2_alphaRho1[celli]   *rAlphaRho1 + (Dinv_T2_alphaRhoU1[celli]    & rAlphaRhoU1) + Dinv_T2_alphaRhoE1[celli]   *rEpsilon1 + Dinv_T2_alphaRho2[celli]   *rAlphaRho2 + (Dinv_T2_alphaRhoU2[celli]    & rAlphaRhoU2) + Dinv_T2_alphaRhoE2[celli]   *rEpsilon2;*/

        //Eigen::Vector<double,10> dState = D_inv[celli].solve((Eigen::Vector<double,10>() << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2).finished());
        
        Eigen::Vector<long double,10> b;
        b << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2;

        Eigen::Vector<long double,10> dState = D_inv[celli].solve(b);

        dp[celli]     = dpStar[celli]     + dState(0);
        dalpha[celli] = dalphaStar[celli] + dState(1);
        dU1[celli]    = dU1Star[celli]    + vector(dState(2), dState(3), dState(4));
        dU2[celli]    = dU2Star[celli]    + vector(dState(5), dState(6), dState(7));
        dT1[celli]    = dT1Star[celli]    + dState(8);
        dT2[celli]    = dT2Star[celli]    + dState(9);
    }
}
