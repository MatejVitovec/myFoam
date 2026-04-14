{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.cells(), celli)
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
                label own = owner[facei];
                label nei = neighbour[facei];

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

                    rAlphaRho1  += 0.5*(                                                   alpha1[own]*unMLambda1*drhodp1[own]*dpStar[own]                              - unMLambda1*rho1[own]*dalphaStar[own]                                    + ((alpha1[own]*rho1[own]*Sf[facei]) & dU1Star[own])                                + alpha1[own]*(unMLambda1*drhodT1[own])*dT1Star[own]);
                    rAlphaRho2  += 0.5*(                                                   alpha2[own]*unMLambda2*drhodp2[own]*dpStar[own]                              + unMLambda2*rho2[own]*dalphaStar[own]                                    + ((alpha2[own]*rho2[own]*Sf[facei]) & dU2Star[own])                                + alpha2[own]*(unMLambda2*drhodT2[own])*dT2Star[own]);

                    rAlphaRhoU1 += 0.5*(                             alpha1[own]*(unMLambda1*drhodp1[own]*U1[own] + Sf[facei])*dpStar[own] - (unMLambda1*rho1[own]*U1[own] + p[own]*Sf[facei])*dalphaStar[own] + (alpha1[own]*rho1[own]*(unMLambda1*tensor::I + (U1[own] * Sf[facei])) & dU1Star[own])                        + alpha1[own]*(unMLambda1*drhodT1[own]*U1[own])*dT1Star[own]);
                    rAlphaRhoU2 += 0.5*(                             alpha2[own]*(unMLambda2*drhodp2[own]*U2[own] + Sf[facei])*dpStar[own] + (unMLambda2*rho2[own]*U2[own] + p[own]*Sf[facei])*dalphaStar[own] + (alpha2[own]*rho2[own]*(unMLambda2*tensor::I + (U2[own] * Sf[facei])) & dU2Star[own])                        + alpha2[own]*(unMLambda2*drhodT2[own]*U2[own])*dT2Star[own]);
                
                    rEpsilon1   += 0.5*(alpha1[own]*(unMLambda1*(drhodp1[own]*H1[own] + drhodT1[own]*T1[own]/rho1[own]) + un1)*dpStar[own]            - (unMLambda1*rho1[own]*E1 + p[own]*un1)*dalphaStar[own]     + ((alpha1[own]*rho1[own]*(unMLambda1*U1[own] + H1[own]*Sf[facei])) & dU1Star[own]) + alpha1[own]*(unMLambda1*(drhodT1[own]*H1[own] + rho1[own]*cp1[own]))*dT1Star[own]);
                    rEpsilon2   += 0.5*(alpha2[own]*(unMLambda2*(drhodp2[own]*H2[own] + drhodT2[own]*T2[own]/rho2[own]) + un2)*dpStar[own]            + (unMLambda2*rho2[own]*E2 + p[own]*un2)*dalphaStar[own]     + ((alpha2[own]*rho2[own]*(unMLambda2*U2[own] + H2[own]*Sf[facei])) & dU2Star[own]) + alpha2[own]*(unMLambda2*(drhodT2[own]*H2[own] + rho2[own]*cp2[own]))*dT2Star[own]);
                }

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

                    rAlphaRho1  += 0.5*(                                                   alpha1[nei]*unMLambda1*drhodp1[nei]*dpStar[nei]                              - unMLambda1*rho1[nei]*dalphaStar[nei]                                    + ((alpha1[nei]*rho1[nei]*Sf[facei]) & dU1Star[nei])                                + alpha1[nei]*(unMLambda1*drhodT1[nei])*dT1Star[nei]);
                    rAlphaRho2  += 0.5*(                                                   alpha2[nei]*unMLambda2*drhodp2[nei]*dpStar[nei]                              + unMLambda2*rho2[nei]*dalphaStar[nei]                                    + ((alpha2[nei]*rho2[nei]*Sf[facei]) & dU2Star[nei])                                + alpha2[nei]*(unMLambda2*drhodT2[nei])*dT2Star[nei]);

                    rAlphaRhoU1 += 0.5*(                             alpha1[nei]*(unMLambda1*drhodp1[nei]*U1[nei] + Sf[facei])*dpStar[nei] - (unMLambda1*rho1[nei]*U1[nei] + p[nei]*Sf[facei])*dalphaStar[nei] + (alpha1[nei]*rho1[nei]*(unMLambda1*tensor::I + (U1[nei] * Sf[facei])) & dU1Star[nei])                        + alpha1[nei]*(unMLambda1*drhodT1[nei]*U1[nei])*dT1Star[nei]);
                    rAlphaRhoU2 += 0.5*(                             alpha2[nei]*(unMLambda2*drhodp2[nei]*U2[nei] + Sf[facei])*dpStar[nei] + (unMLambda2*rho2[nei]*U2[nei] + p[nei]*Sf[facei])*dalphaStar[nei] + (alpha2[nei]*rho2[nei]*(unMLambda2*tensor::I + (U2[nei] * Sf[facei])) & dU2Star[nei])                        + alpha2[nei]*(unMLambda2*drhodT2[nei]*U2[nei])*dT2Star[nei]);
                
                    rEpsilon1   += 0.5*(alpha1[nei]*(unMLambda1*(drhodp1[nei]*H1[nei] + drhodT1[nei]*T1[nei]/rho1[nei]) + un1)*dpStar[nei]            - (unMLambda1*rho1[nei]*E1 + p[nei]*un1)*dalphaStar[nei]     + ((alpha1[nei]*rho1[nei]*(unMLambda1*U1[nei] + H1[nei]*Sf[facei])) & dU1Star[nei]) + alpha1[nei]*(unMLambda1*(drhodT1[nei]*H1[nei] + rho1[nei]*cp1[nei]))*dT1Star[nei]);
                    rEpsilon2   += 0.5*(alpha2[nei]*(unMLambda2*(drhodp2[nei]*H2[nei] + drhodT2[nei]*T2[nei]/rho2[nei]) + un2)*dpStar[nei]            + (unMLambda2*rho2[nei]*E2 + p[nei]*un2)*dalphaStar[nei]     + ((alpha2[nei]*rho2[nei]*(unMLambda2*U2[nei] + H2[nei]*Sf[facei])) & dU2Star[nei]) + alpha2[nei]*(unMLambda2*(drhodT2[nei]*H2[nei] + rho2[nei]*cp2[nei]))*dT2Star[nei]);
                }
            }
            else
            {
                label patchID = mesh.boundaryMesh().whichPatch(facei);
                label localFaceID = facei - mesh.boundaryMesh()[patchID].start();

                //if(!isA<emptyPolyPatch>(mesh.boundaryMesh()[patchID]) && !mesh.boundaryMesh()[patchID].coupled())
                if(isA<wallPolyPatch>(mesh.boundaryMesh()[patchID]))
                {
                    const vector pSf = mesh.Sf().boundaryField()[patchID][localFaceID];
                    const scalar pMagSf = mag(pSf);

                    const vector Uw1 = (tensor::I - 2.0*(pSf * pSf)) & U1[celli];
                    const vector Uw2 = (tensor::I - 2.0*(pSf * pSf)) & U2[celli];

                    //const scalar ac1 = lusgsOmega*(max(max(mag(U1[own] & pSf) + (a1[own] * pMagSf), mag(U2[own] & pSf) + (a2[own]*pMagSf)),
                    //                                   max(mag(U1[celli] & pSf) + (a1[nei] * pMagSf), mag(U2[nei] & pSf) + (a2[nei]*pMagSf))));
                    const scalar ac1 = lusgsOmega*(max(mag(Uw1 & pSf) + (a1[celli] * pMagSf), mag(Uw2 & pSf) + (a2[celli]*pMagSf)));
                    const scalar ac2 = ac1;
                    //const scalar ac1 = lusgsOmega*(mag(U1[celli] & pSf) + a1[celli]*pMagSf);
                    //const scalar ac2 = lusgsOmega*(mag(U2[celli] & pSf) + a2[celli]*pMagSf);

                    //const scalar p_1 = pStar.boundaryField()[patchID][localFaceID];

                    //fluid.blendVanishingFluid(alpha_1, U1_1, U2_1, T1_1, T2_1);
                    //fluid.boundAlpha(alpha_1);

                    const scalar Ew1   = e1[celli] + 0.5*magSqr(Uw1);
                    const scalar Ew2   = e2[celli] + 0.5*magSqr(Uw2);

                    const scalar Hw1 = Ew1 + p[celli]/rho1[celli];
                    const scalar Hw2 = Ew2 + p[celli]/rho2[celli];

                    const scalar un1 = - (U1[celli] & pSf);
                    const scalar un2 = - (U2[celli] & pSf);
                    const scalar unMLambda1 = ac1 - un1;
                    const scalar unMLambda2 = ac2 - un2;

                    const vector dUw1Star = (tensor::I - 2.0*(pSf * pSf)) & dU1Star[celli];
                    const vector dUw2Star = (tensor::I - 2.0*(pSf * pSf)) & dU2Star[celli];

                    rAlphaRho1  += 0.5*(                                                     alpha1[celli]*unMLambda1*drhodp1[celli]*dpStar[celli]                            - unMLambda1*rho1[celli]*dalphaStar[celli]                                + ((alpha1[celli]*rho1[celli]*pSf) & dUw1Star)                                + alpha1[celli]*(unMLambda1*drhodT1[celli])*dT1Star[celli]);
                    rAlphaRho2  += 0.5*(                                                     alpha2[celli]*unMLambda2*drhodp2[celli]*dpStar[celli]                            + unMLambda2*rho2[celli]*dalphaStar[celli]                                + ((alpha2[celli]*rho2[celli]*pSf) & dUw2Star)                                + alpha2[celli]*(unMLambda2*drhodT2[celli])*dT2Star[celli]);

                    rAlphaRhoU1 += 0.5*(                                   alpha1[celli]*(unMLambda1*drhodp1[celli]*Uw1 + pSf)*dpStar[celli] - (unMLambda1*rho1[celli]*Uw1 + p[celli]*pSf)*dalphaStar[celli] + (alpha1[celli]*rho1[celli]*(unMLambda1*tensor::I + (Uw1 * pSf)) & dUw1Star)                            + alpha1[celli]*(unMLambda1*drhodT1[celli]*Uw1)*dT1Star[celli]);
                    rAlphaRhoU2 += 0.5*(                                   alpha2[celli]*(unMLambda2*drhodp2[celli]*Uw2 + pSf)*dpStar[celli] + (unMLambda2*rho2[celli]*Uw2 + p[celli]*pSf)*dalphaStar[celli] + (alpha2[celli]*rho2[celli]*(unMLambda2*tensor::I + (Uw2 * pSf)) & dUw2Star)                            + alpha2[celli]*(unMLambda2*drhodT2[celli]*Uw2)*dT2Star[celli]);
                
                    rEpsilon1   += 0.5*(alpha1[celli]*(unMLambda1*(drhodp1[celli]*Hw1 + drhodT1[celli]*T1[celli]/rho1[celli]) + un1)*dpStar[celli]       - (unMLambda1*rho1[celli]*Ew1 + p[celli]*un1)*dalphaStar[celli]         + ((alpha1[celli]*rho1[celli]*(unMLambda1*Uw1 + Hw1*pSf)) & dUw1Star) + alpha1[celli]*(unMLambda1*(drhodT1[celli]*Hw1 + rho1[celli]*cp1[celli]))*dT1Star[celli]);
                    rEpsilon2   += 0.5*(alpha2[celli]*(unMLambda2*(drhodp2[celli]*Hw2 + drhodT2[celli]*T2[celli]/rho2[celli]) + un2)*dpStar[celli]       + (unMLambda2*rho2[celli]*Ew2 + p[celli]*un2)*dalphaStar[celli]         + ((alpha2[celli]*rho2[celli]*(unMLambda2*Uw2 + Hw2*pSf)) & dUw2Star) + alpha2[celli]*(unMLambda2*(drhodT2[celli]*Hw2 + rho2[celli]*cp2[celli]))*dT2Star[celli]);
                }
            }
        }
    
        Eigen::Vector<double,10> b;
        b << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2;

        //Eigen::Vector<double,10> starState = D_inv[celli].solve((Eigen::Vector<double,10>() << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2).finished());
        Eigen::Vector<double,10> starState = D_inv[celli].solve(b);

        dpStar[celli]     = starState(0);
        dalphaStar[celli] = starState(1);
        dU1Star[celli]    = vector(starState(2), starState(3), starState(4));
        dU2Star[celli]    = vector(starState(5), starState(6), starState(7));
        dT1Star[celli]    = starState(8);
        dT2Star[celli]    = starState(9);
    }
}
