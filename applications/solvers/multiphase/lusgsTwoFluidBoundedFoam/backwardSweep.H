{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), celli)
    {
        const labelList& cellFaces = mesh.cells()[celli];

        scalar rAlphaRho1  = (mesh.V()[celli]/dt.value())*rezAlphaRho1[celli];
        scalar rAlphaRho2  = (mesh.V()[celli]/dt.value())*rezAlphaRho2[celli];
        vector rAlphaRhoU1 = (mesh.V()[celli]/dt.value())*rezAlphaRhoU1[celli];
        vector rAlphaRhoU2 = (mesh.V()[celli]/dt.value())*rezAlphaRhoU2[celli];
        scalar rEpsilon1   = (mesh.V()[celli]/dt.value())*rezEpsilon1[celli];
        scalar rEpsilon2   = (mesh.V()[celli]/dt.value())*rezEpsilon2[celli];
            
        forAll(cellFaces,i)
        {
            label facei = cellFaces[i];

            if (mesh.isInternalFace(facei))
            {
                label nei = neighbour[facei];
                label own = owner[facei];

                if (nei != celli)
                {
                    //const scalar ac1 = lusgsOmega*(max(mag(U1[nei] & Sf[facei]) + (a1[nei] * magSf[facei]), mag(U2[nei] & Sf[facei]) + (a2[nei]*magSf[facei])));
                    //const scalar ac2 = ac1;
                    const scalar ac1 = lusgsOmega*(mag(U1[nei] & Sf[facei]) + a1[nei]*magSf[facei]);
                    const scalar ac2 = lusgsOmega*(mag(U2[nei] & Sf[facei]) + a2[nei]*magSf[facei]);

                    const scalar p_1 = p[nei]     + dp[nei];
                    //scalar alpha_1   = alpha[nei] + dalpha[nei];
                    scalar psi_1     = psi[nei]   + dpsi[nei];
                    vector U1_1      = U1[nei]    + dU1[nei];
                    vector U2_1      = U2[nei]    + dU2[nei];
                    scalar T1_1      = T1[nei]    + dT1[nei];
                    scalar T2_1      = T2[nei]    + dT2[nei];

                    //fluid.blendVanishingFluid(alpha_1, U1_1, U2_1, T1_1, T2_1);
                    //fluid.boundAlpha(alpha_1);
                    psi_1 = min(psi_1, 50.0);
                    psi_1 = max(psi_1, -50.0);

                    const scalar alpha2_1 = 1.0/(1.0 + Foam::exp(-psi_1));
                    const scalar alpha1_1 = 1.0 - alpha2_1;
                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);

                    const scalar E1   = e1[nei]                 + 0.5*magSqr(U1[nei]);
                    const scalar E2   = e2[nei]                 + 0.5*magSqr(U2[nei]);
                    const scalar E1_1 = gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1);
                    const scalar E2_1 = gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1);

                    const scalar alphaRhoPhi1   = alpha1[nei]*rho1[nei]*(U1[nei] & Sf[facei]);
                    const scalar alphaRhoPhi2   = alpha2[nei]*rho2[nei]*(U2[nei] & Sf[facei]);
                    const scalar alphaRhoPhi1_1 = alpha1_1*rho1_1*(U1_1 & Sf[facei]);
                    const scalar alphaRhoPhi2_1 = alpha2_1*rho2_1*(U2_1 & Sf[facei]);

                    rAlphaRho1  += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1      - alpha1[nei]*rho1[nei])         - 0.5*(alphaRhoPhi1_1                                 - alphaRhoPhi1));
                    rAlphaRho2  += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1      - alpha2[nei]*rho2[nei])         - 0.5*(alphaRhoPhi2_1                                 - alphaRhoPhi2)); 
                    rAlphaRhoU1 += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1*U1_1 - alpha1[nei]*rho1[nei]*U1[nei]) - 0.5*((alphaRhoPhi1_1*U1_1 + alpha1_1*p_1*Sf[facei]) - (alphaRhoPhi1*U1[nei] + alpha1[nei]*p[nei]*Sf[facei])));
                    rAlphaRhoU2 += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1*U2_1 - alpha2[nei]*rho2[nei]*U2[nei]) - 0.5*((alphaRhoPhi2_1*U2_1 + alpha2_1*p_1*Sf[facei]) - (alphaRhoPhi2*U2[nei] + alpha2[nei]*p[nei]*Sf[facei])));
                    rEpsilon1   += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1*E1_1 - alpha1[nei]*rho1[nei]*E1)      - 0.5*(alphaRhoPhi1_1*(E1_1 + p_1/rho1_1)             - alphaRhoPhi1*(E1 + p[nei]/rho1[nei])));
                    rEpsilon2   += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1*E2_1 - alpha2[nei]*rho2[nei]*E2)      - 0.5*(alphaRhoPhi2_1*(E2_1 + p_1/rho2_1)             - alphaRhoPhi2*(E2 + p[nei]/rho2[nei])));
                }

                if (own != celli)
                {
                    //const scalar ac1 = lusgsOmega*(max(mag(U1[own] & Sf[facei]) + (a1[own] * magSf[facei]), mag(U2[own] & Sf[facei]) + (a2[own]*magSf[facei])));
                    //const scalar ac2 = ac1;
                    const scalar ac1 = lusgsOmega*(mag(U1[own] & Sf[facei]) + a1[own]*magSf[facei]);
                    const scalar ac2 = lusgsOmega*(mag(U2[own] & Sf[facei]) + a2[own]*magSf[facei]);
                    
                    const scalar p_1 = p[own]     + dpStar[own];
                    //scalar alpha_1   = alpha[own] + dalphaStar[own];
                    scalar psi_1     = psi[own]   + dpsiStar[own];
                    vector U1_1      = U1[own]    + dU1Star[own];
                    vector U2_1      = U2[own]    + dU2Star[own];
                    scalar T1_1      = T1[own]    + dT1Star[own];
                    scalar T2_1      = T2[own]    + dT2Star[own];

                    //fluid.blendVanishingFluid(alpha_1, U1_1, U2_1, T1_1, T2_1);
                    //fluid.boundAlpha(alpha_1);
                    psi_1 = min(psi_1, 50.0);
                    psi_1 = max(psi_1, -50.0);

                    const scalar alpha2_1 = 1.0/(1.0 + Foam::exp(-psi_1));
                    const scalar alpha1_1 = 1.0 - alpha2_1;
                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);

                    const scalar E1   = e1[own]                 + 0.5*magSqr(U1[own]);
                    const scalar E2   = e2[own]                 + 0.5*magSqr(U2[own]);
                    const scalar E1_1 = gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1);
                    const scalar E2_1 = gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1);

                    const scalar alphaRhoPhi1   = alpha1[own]*rho1[own]*(U1[own] & Sf[facei]);
                    const scalar alphaRhoPhi2   = alpha2[own]*rho2[own]*(U2[own] & Sf[facei]);
                    const scalar alphaRhoPhi1_1 = alpha1_1*rho1_1*(U1_1 & Sf[facei]);
                    const scalar alphaRhoPhi2_1 = alpha2_1*rho2_1*(U2_1 & Sf[facei]);

                    rAlphaRho1  += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1      - alpha1[own]*rho1[own])         + 0.5*(alphaRhoPhi1_1                                 - alphaRhoPhi1));
                    rAlphaRho2  += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1      - alpha2[own]*rho2[own])         + 0.5*(alphaRhoPhi2_1                                 - alphaRhoPhi2)); 
                    rAlphaRhoU1 += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1*U1_1 - alpha1[own]*rho1[own]*U1[own]) + 0.5*((alphaRhoPhi1_1*U1_1 + alpha1_1*p_1*Sf[facei]) - (alphaRhoPhi1*U1[own] + alpha1[own]*p[own]*Sf[facei])));
                    rAlphaRhoU2 += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1*U2_1 - alpha2[own]*rho2[own]*U2[own]) + 0.5*((alphaRhoPhi2_1*U2_1 + alpha2_1*p_1*Sf[facei]) - (alphaRhoPhi2*U2[own] + alpha2[own]*p[own]*Sf[facei])));
                    rEpsilon1   += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1*E1_1 - alpha1[own]*rho1[own]*E1)      + 0.5*(alphaRhoPhi1_1*(E1_1 + p_1/rho1_1)             - alphaRhoPhi1*(E1 + p[own]/rho1[own])));
                    rEpsilon2   += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1*E2_1 - alpha2[own]*rho2[own]*E2)      + 0.5*(alphaRhoPhi2_1*(E2_1 + p_1/rho2_1)             - alphaRhoPhi2*(E2 + p[own]/rho2[own])));
                }
            }
            /*else
            {
                label patchID = mesh.boundaryMesh().whichPatch(facei);
                label localFaceID = facei - mesh.boundaryMesh()[patchID].start();

                if(!isA<emptyPolyPatch>(mesh.boundaryMesh()[patchID]) && !mesh.boundaryMesh()[patchID].coupled())
                {
                    const vector pSf = mesh.Sf().boundaryField()[patchID][localFaceID];
                    const scalar pMagSf = mag(pSf);

                    const scalar p_1 = pStar.boundaryField()[patchID][localFaceID];
                    scalar alpha_1   = alphaStar.boundaryField()[patchID][localFaceID];
                    vector U1_1      = U1Star.boundaryField()[patchID][localFaceID];
                    vector U2_1      = U2Star.boundaryField()[patchID][localFaceID];
                    scalar T1_1      = T1Star.boundaryField()[patchID][localFaceID];
                    scalar T2_1      = T2Star.boundaryField()[patchID][localFaceID];

                    const scalar p_0        = p.boundaryField()[patchID][localFaceID];
                    const scalar alpha_0    = alpha.boundaryField()[patchID][localFaceID];
                    const vector U1_0       = U1.boundaryField()[patchID][localFaceID];
                    const vector U2_0       = U2.boundaryField()[patchID][localFaceID];
                    const scalar T1_0       = T1.boundaryField()[patchID][localFaceID];
                    const scalar T2_0       = T2.boundaryField()[patchID][localFaceID];
                    const scalar a1_0       = a1.boundaryField()[patchID][localFaceID];
                    const scalar a2_0       = a2.boundaryField()[patchID][localFaceID];

                    //const scalar ac1 = lusgsOmega*(max(mag(U1_0 & pSf) + (a1_0*pMagSf), mag(U2_0 & pSf) + (a2_0*pMagSf)));
                    //const scalar ac2 = ac1;
                    const scalar ac1 = lusgsOmega*(mag(U1_0 & pSf) + a1_0*pMagSf);
                    const scalar ac2 = lusgsOmega*(mag(U2_0 & pSf) + a2_0*pMagSf);

                    //fluid.blendVanishingFluid(alpha_1, U1_1, U2_1, T1_1, T2_1);
                    //fluid.boundAlpha(alpha_1);

                    const scalar alpha1_0 = 1.0 - alpha_0;
                    const scalar alpha2_0 = alpha_0;
                    const scalar alpha1_1 = 1.0 - alpha_1;
                    const scalar alpha2_1 = alpha_1;

                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);
                    const scalar rho1_0 = gasProps1.rho(p_0, T1_0);
                    const scalar rho2_0 = gasProps2.rho(p_0, T2_0);

                    const scalar E1_0 = gasProps1.Es(p_0, T1_0) + 0.5*magSqr(U1_0);
                    const scalar E2_0 = gasProps2.Es(p_0, T2_0) + 0.5*magSqr(U2_0);
                    const scalar E1_1 = gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1);
                    const scalar E2_1 = gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1);

                    const scalar alphaRhoPhi1_0 = alpha1_0*rho1_0*(U1_0 & pSf);
                    const scalar alphaRhoPhi2_0 = alpha2_0*rho2_0*(U2_0 & pSf);
                    const scalar alphaRhoPhi1_1 = alpha1_1*rho1_1*(U1_1 & pSf);
                    const scalar alphaRhoPhi2_1 = alpha2_1*rho2_1*(U2_1 & pSf);

                    rAlphaRho1  += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1      - alpha1_0*rho1_0)         + 0.5*(alphaRhoPhi1_1                           - alphaRhoPhi1_0));
                    rAlphaRho2  += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1      - alpha2_0*rho2_0)         + 0.5*(alphaRhoPhi2_1                           - alphaRhoPhi2_0)); 
                    rAlphaRhoU1 += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1*U1_1 - alpha1_0*rho1_0*U1_0)    + 0.5*((alphaRhoPhi1_1*U1_1 + alpha1_1*p_1*pSf) - (alphaRhoPhi1_0*U1_0 + alpha1_0*p_0*pSf)));
                    rAlphaRhoU2 += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1*U2_1 - alpha2_0*rho2_0*U2_0)    + 0.5*((alphaRhoPhi2_1*U2_1 + alpha2_1*p_1*pSf) - (alphaRhoPhi2_0*U2_0 + alpha2_0*p_0*pSf)));
                    rEpsilon1   += lusgsRelax*(0.5*ac1*(alpha1_1*rho1_1*E1_1 - alpha1_0*rho1_0*E1_0)    + 0.5*(alphaRhoPhi1_1*(E1_1 + p_1/rho1_1)       - alphaRhoPhi1_0*(E1_0 + p_0/rho1_0)));
                    rEpsilon2   += lusgsRelax*(0.5*ac2*(alpha2_1*rho2_1*E2_1 - alpha2_0*rho2_0*E2_0)    + 0.5*(alphaRhoPhi2_1*(E2_1 + p_1/rho2_1)       - alphaRhoPhi2_0*(E2_0 + p_0/rho2_0)));
                }
            }*/
        }
        
        Eigen::Vector<double,10> b;
        b << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2;

        //Eigen::Vector<double,10> dState = D_inv[celli].solve((Eigen::Vector<double,10>() << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2).finished()); 
        Eigen::Vector<double,10> dState = D_inv[celli].solve(b);

        dp[celli]     = dState(0);
        //dalpha[celli] = dState(1);
        dpsi[celli] = dState(1);
        dU1[celli]    = vector(dState(2), dState(3), dState(4));
        dU2[celli]    = vector(dState(5), dState(6), dState(7));
        dT1[celli]    = dState(8);
        dT2[celli]    = dState(9);
    }
}
