{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.cells(), celli)
    {
        const labelList& cellFaces = mesh.cells()[celli];

        scalar rRho  = mesh.V()[celli]/dt.value()*resRho[celli];
        vector rRhoU = mesh.V()[celli]/dt.value()*resRhoU[celli];
        scalar rRhoE = mesh.V()[celli]/dt.value()*resRhoE[celli];
        scalar rRhow = mesh.V()[celli]/dt.value()*resRhow[celli];

        forAll(cellFaces, i)
        {
            label faceI = cellFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // we use the mesh ordering property: own < nei
                if (own != celli)  // then data[own] belongs to lower triangle
                {
                    scalar ac = lusgsOmega*(mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);

                    scalar p1 = p[own] + dpStar[own];
                    vector U1 = U[own] + dUStar[own];
                    scalar T1 = T[own] + dTStar[own];
                    scalar w1 = w[own] + dwStar[own];
                    
                    scalar rho1 = 1.0/(w1/gasPropsl.rho(p1, T1) + (1.0 - w1)/gasPropsv.rho(p1, T1));
                    scalar rhoE1 = rho1*(w1*gasPropsl.Es(p1, T1) + (1.0 - w1)*gasPropsv.Es(p1, T1) + 0.5*magSqr(U1));

                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];

                    rRho  += 0.5*ac*(rho1    - rho[own])        + 0.5*(rho1*phi1 - rho[own]*phi0); 
                    rRhoU += 0.5*ac*(rho1*U1 - rho[own]*U[own]) + 0.5*((rho1*U1*phi1 - rho[own]*U[own]*phi0) + (p1 - p[own])*Sf[faceI]);
                    rRhoE += 0.5*ac*(rhoE1   - rho[own]*E[own]) + 0.5*(rhoE1*phi1 + p1*phi1 - (rho[own]*E[own]*phi0 + p[own]*phi0));
                    rRhow += 0.5*ac*(rho1*w1 - rho[own]*w[own]) + 0.5*(rho1*w1*phi1 - rho[own]*w[own]*phi0); 

                    /*scalar ac = lusgsOmega*(mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);

                    scalar un = U[own] & Sf[faceI];
                    scalar dun = dUStar[own] & Sf[faceI];
                    scalar unMac = un + ac;
                    
                    rRho  += 0.5*(unMac*drhodp[own]*dpStar[own]                                                           + rho[own]*dun                                                + unMac*drhodT[own]*dTStar[own]);
                    rRhoU += 0.5*(unMac*drhodp[own]*dpStar[own]*U[own]           + Sf[faceI]*dpStar[own]                  + rho[own]*dun*U[own] + unMac*rho[own]*dUStar[own]            + unMac*drhodT[own]*dTStar[own]*U[own]);
                    rRhoE += 0.5*(unMac*(drhodp[own]*H[own] + (T[own]/rho[own])*drhodT[own])*dpStar[own] + un*dpStar[own] + rho[own]*dun*H[own] + unMac*rho[own]*(U[own] & dUStar[own]) + unMac*(drhodT[own]*H[own] + rho[own]*cp[own])*dTStar[own]);*/
                }

                if (nei != celli)  // then data[own] belongs to lower triangle
                {
                    scalar ac = lusgsOmega*(mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar p1 = p[nei] + dp[nei];
                    vector U1 = U[nei] + dU[nei];
                    scalar T1 = T[nei] + dT[nei];
                    scalar w1 = w[nei] + dwStar[nei];
                    
                    scalar rho1 = 1.0/(w1/gasPropsl.rho(p1, T1) + (1.0 - w1)/gasPropsv.rho(p1, T1));
                    scalar rhoE1 = rho1*(w1*gasPropsl.Es(p1, T1) + (1.0 - w1)*gasPropsv.Es(p1, T1) + 0.5*magSqr(U1));

                    scalar phi0 = U[nei] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];

                    rRho  += 0.5*ac*(rho1    - rho[nei])        - 0.5*(rho1*phi1 - rho[nei]*phi0); 
                    rRhoU += 0.5*ac*(rho1*U1 - rho[nei]*U[nei]) - 0.5*((rho1*U1*phi1 - rho[nei]*U[nei]*phi0) + (p1 - p[nei])*Sf[faceI]);
                    rRhoE += 0.5*ac*(rhoE1   - rho[nei]*E[nei]) - 0.5*(rhoE1*phi1 + p1*phi1 - (rho[nei]*E[nei]*phi0 + p[nei]*phi0));
                    rRhow += 0.5*ac*(rho1*w1 - rho[nei]*w[nei]) - 0.5*(rho1*w1*phi1 - rho[nei]*w[nei]*phi0); 
                    
                    /*scalar ac = lusgsOmega*(mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar un = U[nei] & Sf[faceI];
                    scalar dun = dU[nei] & Sf[faceI];
                    scalar unMac = un - ac;
                    
                    rRho  -= 0.5*(unMac*drhodp[nei]*dp[nei]                                                       + rho[nei]*dun                                            + unMac*drhodT[nei]*dT[nei]);
                    rRhoU -= 0.5*(unMac*drhodp[nei]*dp[nei]*U[nei]               + Sf[faceI]*dp[nei]              + rho[nei]*dun*U[nei] + unMac*rho[nei]*dU[nei]            + unMac*drhodT[nei]*dT[nei]*U[nei]);
                    rRhoE -= 0.5*(unMac*(drhodp[nei]*H[nei] + (T[nei]/rho[nei])*drhodT[nei])*dp[nei] + un*dp[nei] + rho[nei]*dun*H[nei] + unMac*rho[nei]*(U[nei] & dU[nei]) + unMac*(drhodT[nei]*H[nei] + rho[nei]*cp[nei])*dT[nei]);*/
                }
            }
        }

        dpStar[celli] = ((invA[celli]/invD[celli])*rRho + (((drhodT[celli]/invD[celli])*U[celli]) & rRhoU) + (-drhodT[celli]/invD[celli])*rRhoE + (invC[celli]/invD[celli])*rRhow)/D[celli];
        dUStar[celli] = ((-U[celli]/rho[celli])   *rRho + (1/rho[celli])*rRhoU)/D[celli];
        dpStar[celli] = ((invB[celli]/invD[celli])*rRho + (((-drhodp[celli]/invD[celli])*U[celli]) & rRhoU) + (drhodp[celli]/invD[celli])*rRhoE + (invE[celli]/invD[celli])*rRhow)/D[celli];
        dwStar[celli] = ((-w[celli]/rho[celli])   *rRho + (1.0/rho[celli])*rRhow)/D[celli];
    }
}
