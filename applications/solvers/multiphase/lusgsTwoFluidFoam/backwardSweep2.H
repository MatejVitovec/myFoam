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
                    //const scalar ac1 = lusgsOmega*(max(max(mag(U1[own] & Sf[facei]) + (a1[own] * magSf[facei]), mag(U2[own] & Sf[facei]) + (a2[own]*magSf[facei])),
                    //                                   max(mag(U1[nei] & Sf[facei]) + (a1[nei] * magSf[facei]), mag(U2[nei] & Sf[facei]) + (a2[nei]*magSf[facei]))));
                    const scalar ac1 = lusgsOmega*(max(mag(U1[nei] & Sf[facei]) + (a1[nei] * magSf[facei]), mag(U2[nei] & Sf[facei]) + (a2[nei]*magSf[facei])));
                    const scalar ac2 = ac1;
                    //const scalar ac1 = lusgsOmega*(mag(U1[nei] & Sf[facei]) + a1[nei]*magSf[facei]);
                    //const scalar ac2 = lusgsOmega*(mag(U2[nei] & Sf[facei]) + a2[nei]*magSf[facei]);

                    const scalar E1   = e1[nei] + 0.5*magSqr(U1[nei]);
                    const scalar E2   = e2[nei] + 0.5*magSqr(U2[nei]);

                    const scalar un1 = U1[nei] & Sf[facei];
                    const scalar un2 = U2[nei] & Sf[facei];
                    const scalar unMLambda1 = ac1 - un1;
                    const scalar unMLambda2 = ac2 - un2;

                    rAlphaRho1  += 0.5*(                                                   alpha1[nei]*unMLambda1*drhodp1[nei]*dpStar[nei]                              - unMLambda1*rho1[nei]*dalphaStar[nei]                                    + (alpha1[nei]*rho1[nei]*Sf[facei]) & dU1Star[nei]                                + alpha1[nei]*(unMLambda1*drhodT1[nei])*dT1Star[nei]);
                    rAlphaRho2  += 0.5*(                                                   alpha2[nei]*unMLambda2*drhodp2[nei]*dpStar[nei]                              + unMLambda2*rho2[nei]*dalphaStar[nei]                                    + (alpha2[nei]*rho2[nei]*Sf[facei]) & dU2Star[nei]                                + alpha2[nei]*(unMLambda2*drhodT2[nei])*dT2Star[nei]);

                    rAlphaRhoU1 += 0.5*(                             alpha1[nei]*(unMLambda1*drhodp1[nei]*U1[nei] + Sf[facei])*dpStar[nei] - (unMLambda1*rho1[nei]*U1[nei] + p[nei]*Sf[facei])*dalphaStar[nei] + alpha1[nei]*rho1[nei]*(unMLambda1*tensor::I + (U1[nei] * Sf[facei])) & dU1Star[nei]                        + alpha1[nei]*(unMLambda1*drhodT1[nei]*U1[nei])*dT1Star[nei]);
                    rAlphaRhoU2 += 0.5*(                             alpha2[nei]*(unMLambda2*drhodp2[nei]*U1[nei] + Sf[facei])*dpStar[nei] + (unMLambda2*rho2[nei]*U2[nei] + p[nei]*Sf[facei])*dalphaStar[nei] + alpha2[nei]*rho2[nei]*(unMLambda2*tensor::I + (U2[nei] * Sf[facei])) & dU2Star[nei]                        + alpha2[nei]*(unMLambda2*drhodT2[nei]*U2[nei])*dT2Star[nei]);
                
                    rEpsilon1   += 0.5*(alpha1[nei]*(unMLambda1*(drhodp1[nei]*H1[nei] + drhodT1[nei]*T1[nei]/rho1[nei]) + un1)*dpStar[nei]            - (unMLambda1*rho1[nei]*E1 + p[nei]*un1)*dalphaStar[nei]       + alpha1[nei]*rho1[nei]*(unMLambda1*U1[nei] + H1[nei]*Sf[facei]) & dU1Star[nei] + alpha1[nei]*(unMLambda1*(drhodT1[nei]*H1[nei] + rho1[nei]*cp1[nei]))*dT1Star[nei]);
                    rEpsilon2   += 0.5*(alpha2[nei]*(unMLambdaě*(drhodp2[nei]*H2[nei] + drhodT2[nei]*T2[nei]/rho2[nei]) + un2)*dpStar[nei]            + (unMLambda2*rho2[nei]*E2 + p[nei]*un2)*dalphaStar[nei]       + alpha2[nei]*rho2[nei]*(unMLambda2*U2[nei] + H2[nei]*Sf[facei]) & dU2Star[nei] + alpha2[nei]*(unMLambda2*(drhodT2[nei]*H2[nei] + rho2[nei]*cp2[nei]))*dT2Star[nei]);
                }

                if (own != celli)
                {
                    //const scalar ac1 = lusgsOmega*(max(max(mag(U1[own] & Sf[facei]) + (a1[own] * magSf[facei]), mag(U2[own] & Sf[facei]) + (a2[own]*magSf[facei])),
                    //                                   max(mag(U1[nei] & Sf[facei]) + (a1[nei] * magSf[facei]), mag(U2[nei] & Sf[facei]) + (a2[nei]*magSf[facei]))));
                    const scalar ac1 = lusgsOmega*(max(mag(U1[own] & Sf[facei]) + (a1[own] * magSf[facei]), mag(U2[own] & Sf[facei]) + (a2[own]*magSf[facei])));
                    const scalar ac2 = ac1;
                    //const scalar ac1 = lusgsOmega*(mag(U1[own] & Sf[facei]) + a1[own]*magSf[facei]);
                    //const scalar ac2 = lusgsOmega*(mag(U2[own] & Sf[facei]) + a2[own]*magSf[facei]);

                    const scalar E1   = e1[own] + 0.5*magSqr(U1[own]);
                    const scalar E2   = e2[own] + 0.5*magSqr(U2[own]);

                    const scalar un1 = U1[own] & Sf[facei];
                    const scalar un2 = U2[own] & Sf[facei];
                    const scalar unMLambda1 = ac1 + un1;
                    const scalar unMLambda2 = ac2 + un2;

                    rAlphaRho1  += 0.5*(                                                   alpha1[own]*unMLambda1*drhodp1[own]*dpStar[own]                              - unMLambda1*rho1[own]*dalphaStar[own]                                    + (alpha1[own]*rho1[own]*Sf[facei]) & dU1Star[own]                                + alpha1[own]*(unMLambda1*drhodT1[own])*dT1Star[own]);
                    rAlphaRho2  += 0.5*(                                                   alpha2[own]*unMLambda2*drhodp2[own]*dpStar[own]                              + unMLambda2*rho2[own]*dalphaStar[own]                                    + (alpha2[own]*rho2[own]*Sf[facei]) & dU2Star[own]                                + alpha2[own]*(unMLambda2*drhodT2[own])*dT2Star[own]);

                    rAlphaRhoU1 += 0.5*(                             alpha1[own]*(unMLambda1*drhodp1[own]*U1[own] + Sf[facei])*dpStar[own] - (unMLambda1*rho1[own]*U1[own] + p[own]*Sf[facei])*dalphaStar[own] + alpha1[own]*rho1[own]*(unMLambda1*tensor::I + (U1[own] * Sf[facei])) & dU1Star[own]                        + alpha1[own]*(unMLambda1*drhodT1[own]*U1[own])*dT1Star[own]);
                    rAlphaRhoU2 += 0.5*(                             alpha2[own]*(unMLambda2*drhodp2[own]*U1[own] + Sf[facei])*dpStar[own] + (unMLambda2*rho2[own]*U2[own] + p[own]*Sf[facei])*dalphaStar[own] + alpha2[own]*rho2[own]*(unMLambda2*tensor::I + (U2[own] * Sf[facei])) & dU2Star[own]                        + alpha2[own]*(unMLambda2*drhodT2[own]*U2[own])*dT2Star[own]);
                
                    rEpsilon1   += 0.5*(alpha1[own]*(unMLambda1*(drhodp1[own]*H1[own] + drhodT1[own]*T1[own]/rho1[own]) + un1)*dpStar[own]            - (unMLambda1*rho1[own]*E1 + p[own]*un1)*dalphaStar[own]       + alpha1[own]*rho1[own]*(unMLambda1*U1[own] + H1[own]*Sf[facei]) & dU1Star[own] + alpha1[own]*(unMLambda1*(drhodT1[own]*H1[own] + rho1[own]*cp1[own]))*dT1Star[own]);
                    rEpsilon2   += 0.5*(alpha2[own]*(unMLambdaě*(drhodp2[own]*H2[own] + drhodT2[own]*T2[own]/rho2[own]) + un2)*dpStar[own]            + (unMLambda2*rho2[own]*E2 + p[own]*un2)*dalphaStar[own]       + alpha2[own]*rho2[own]*(unMLambda2*U2[own] + H2[own]*Sf[facei]) & dU2Star[own] + alpha2[own]*(unMLambda2*(drhodT2[own]*H2[own] + rho2[own]*cp2[own]))*dT2Star[own]);
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

                    const scalar ac1 = lusgsOmega*(max(mag(U1_0 & pSf) + (a1_0*pMagSf), mag(U2_0 & pSf) + (a2_0*pMagSf)));
                    const scalar ac2 = ac1;
                    //const scalar ac1 = lusgsOmega*(mag(U1_0 & pSf) + a1_0*pMagSf);
                    //const scalar ac2 = lusgsOmega*(mag(U2_0 & pSf) + a2_0*pMagSf);

                    //fluid.blendVanishingFluid(alpha_1, U1_1, U2_1, T1_1, T2_1);
                    //fluid.boundAlpha(alpha_1);

                    const scalar alpha1_0 = 1.0 - alpha_0;
                    const scalar alpha2_0 = alpha_0;
                    const scalar alpha1_1 = 1.0 - alpha_1;
                    const scalar alpha2_1 = alpha_1;

                    const scalar rho1_0 = gasProps1.rho(p_0, T1_0);
                    const scalar rho2_0 = gasProps2.rho(p_0, T2_0);
                    const scalar rho1_1 = gasProps1.rho(p_1, T1_1);
                    const scalar rho2_1 = gasProps2.rho(p_1, T2_1);

                    const scalar E1_0 = gasProps1.Es(p_0, T1_0) + 0.5*magSqr(U1_0);
                    const scalar E2_0 = gasProps2.Es(p_0, T2_0) + 0.5*magSqr(U2_0);
                    const scalar E1_1 = gasProps1.Es(p_1, T1_1) + 0.5*magSqr(U1_1);
                    const scalar E2_1 = gasProps2.Es(p_1, T2_1) + 0.5*magSqr(U2_1);

                    const scalar alphaRhoPhi1_0 = alpha1_0*rho1_0*(U1_0 & pSf);
                    const scalar alphaRhoPhi2_0 = alpha2_0*rho2_0*(U2_0 & pSf);
                    const scalar alphaRhoPhi1_1 = alpha1_1*rho1_1*(U1_1 & pSf);
                    const scalar alphaRhoPhi2_1 = alpha2_1*rho2_1*(U2_1 & pSf);

                    rAlphaRho1  += 0.5*ac1*(alpha1_1*rho1_1      - alpha1_0*rho1_0)         - 0.5*(alphaRhoPhi1_1                           - alphaRhoPhi1_0);
                    rAlphaRho2  += 0.5*ac2*(alpha2_1*rho2_1      - alpha2_0*rho2_0)         - 0.5*(alphaRhoPhi2_1                           - alphaRhoPhi2_0); 
                    rAlphaRhoU1 += 0.5*ac1*(alpha1_1*rho1_1*U1_1 - alpha1_0*rho1_0*U1_0)    - 0.5*((alphaRhoPhi1_1*U1_1 + alpha1_1*p_1*pSf) - (alphaRhoPhi1_0*U1_0 + alpha1_0*p_0*pSf));
                    rAlphaRhoU2 += 0.5*ac2*(alpha2_1*rho2_1*U2_1 - alpha2_0*rho2_0*U2_0)    - 0.5*((alphaRhoPhi2_1*U2_1 + alpha2_1*p_1*pSf) - (alphaRhoPhi2_0*U2_0 + alpha2_0*p_0*pSf));
                    rEpsilon1   += 0.5*ac1*(alpha1_1*rho1_1*E1_1 - alpha1_0*rho1_0*E1_0)    - 0.5*(alphaRhoPhi1_1*(E1_1 + p_1/rho1_1)       - alphaRhoPhi1_0*(E1_0 + p_0/rho1_0));
                    rEpsilon2   += 0.5*ac2*(alpha2_1*rho2_1*E2_1 - alpha2_0*rho2_0*E2_0)    - 0.5*(alphaRhoPhi2_1*(E2_1 + p_1/rho2_1)       - alphaRhoPhi2_0*(E2_0 + p_0/rho2_0));
                }
            }*/
        }
        
        Eigen::Vector<double,10> b;
        b << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2;

        //Eigen::Vector<double,10> dState = D_inv[celli].solve((Eigen::Vector<double,10>() << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2).finished()); 
        Eigen::Vector<double,10> dState = D_inv[celli].solve(b);

        dp[celli]     = dState(0);
        dalpha[celli] = dState(1);
        dU1[celli]    = vector(dState(2), dState(3), dState(4));
        dU2[celli]    = vector(dState(5), dState(6), dState(7));
        dT1[celli]    = dState(8);
        dT2[celli]    = dState(9);
    }
}
