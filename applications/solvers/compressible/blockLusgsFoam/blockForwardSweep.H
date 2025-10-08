{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];

        scalar rRho  = mesh.V()[cellI]/dt.value()*resRho[cellI];
        vector rRhoU = mesh.V()[cellI]/dt.value()*resRhoU[cellI];
        scalar rRhoE = mesh.V()[cellI]/dt.value()*resRhoE[cellI];

        forAll(cellFaces, i)
        {
            label faceI = cellFaces[i];

            if (mesh.isInternalFace(faceI))
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                    /*scalar ac = lusgsOmega*(mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);

                    scalar rhoE = rho[own]*(e[own] + 0.5*magSqr(U[own]));

                    scalar p1 = p[own] + dpStar[own];
                    vector U1 = U[own] + dUStar[own];
                    scalar T1 = T[own] + dTStar[own];
                    
                    scalar rho1 = gasProps.rho(p1, T1);
                    scalar rhoE1 = rho1*(gasProps.Es(p1, T1) + 0.5*magSqr(U1));

                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];

                    rRho  += 0.5*ac*(rho1    - rho[own])        + 0.5*(rho1*phi1 - rho[own]*phi0); 
                    rRhoU += 0.5*ac*(rho1*U1 - rho[own]*U[own]) + 0.5*((rho1*U1*phi1 - rho[own]*U[own]*phi0) + (p1 - p[own])*Sf[faceI]);
                    rRhoE += 0.5*ac*(rhoE1   - rhoE)            + 0.5*(rhoE1*phi1 + p1*phi1 - (rhoE*phi0 + p[own]*phi0));*/

                    scalar ac = lusgsOmega*(mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);

                    scalar un = U[own] & Sf[faceI];
                    scalar dun = dUStar[own] & Sf[faceI];
                    scalar unMac = un + ac;
                    
                    rRho  += 0.5*(unMac*drhodp[own]*dpStar[own]                                                           + rho[own]*dun                                                + unMac*drhodT[own]*dTStar[own]);
                    rRhoU += 0.5*(unMac*drhodp[own]*dpStar[own]*U[own]           + Sf[faceI]*dpStar[own]                  + rho[own]*dun*U[own] + unMac*rho[own]*dUStar[own]            + unMac*drhodT[own]*dTStar[own]*U[own]);
                    rRhoE += 0.5*(unMac*(drhodp[own]*H[own] + (T[own]/rho[own])*drhodT[own])*dpStar[own] + un*dpStar[own] + rho[own]*dun*H[own] + unMac*rho[own]*(U[own] & dUStar[own]) + unMac*(drhodT[own]*H[own] + rho[own]*cp[own])*dTStar[own]);
                }

                if (nei != cellI)  // then data[own] belongs to lower triangle
                {
                    /*scalar ac = lusgsOmega*(mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar rhoE = rho[nei]*(e[nei] + 0.5*magSqr(U[nei]));

                    scalar p1 = p[nei] + dp[nei];
                    vector U1 = U[nei] + dU[nei];
                    scalar T1 = T[nei] + dT[nei];
                    
                    scalar rho1 = gasProps.rho(p1, T1);
                    scalar rhoE1 = rho1*(gasProps.Es(p1, T1) + 0.5*magSqr(U1));

                    scalar phi0 = U[nei] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];

                    rRho  += 0.5*ac*(rho1    - rho[nei])        - 0.5*(rho1*phi1 - rho[nei]*phi0); 
                    rRhoU += 0.5*ac*(rho1*U1 - rho[nei]*U[nei]) - 0.5*((rho1*U1*phi1 - rho[nei]*U[nei]*phi0) + (p1 - p[nei])*Sf[faceI]);
                    rRhoE += 0.5*ac*(rhoE1   - rhoE)            - 0.5*(rhoE1*phi1 + p1*phi1 - (rhoE*phi0 + p[nei]*phi0));*/
                    
                    scalar ac = lusgsOmega*(mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar un = U[nei] & Sf[faceI];
                    scalar dun = dU[nei] & Sf[faceI];
                    scalar unMac = un - ac;
                    
                    rRho  -= 0.5*(unMac*drhodp[nei]*dp[nei]                                                           + rho[nei]*dun                                            + unMac*drhodT[nei]*dT[nei]);
                    rRhoU -= 0.5*(unMac*drhodp[nei]*dp[nei]*U[nei]               + Sf[faceI]*dp[nei]                  + rho[nei]*dun*U[nei] + unMac*rho[nei]*dU[nei]            + unMac*drhodT[nei]*dT[nei]*U[nei]);
                    rRhoE -= 0.5*(unMac*(drhodp[nei]*H[nei] + (T[nei]/rho[nei])*drhodT[nei])*dp[nei] + un*dp[nei]     + rho[nei]*dun*H[nei] + unMac*rho[nei]*(U[nei] & dU[nei]) + unMac*(drhodT[nei]*H[nei] + rho[nei]*cp[nei])*dT[nei]);
                }
            }
        }

        dpStar[cellI] = Dinv_p_rho[cellI]*rRho + (Dinv_p_rhoU[cellI] & rRhoU) + Dinv_p_rhoE[cellI]*rRhoE;
        dUStar[cellI] = Dinv_U_rho[cellI]*rRho + (Dinv_U_rhoU[cellI] & rRhoU) + Dinv_U_rhoE[cellI]*rRhoE;
        dTStar[cellI] = Dinv_T_rho[cellI]*rRho + (Dinv_T_rhoU[cellI] & rRhoU) + Dinv_T_rhoE[cellI]*rRhoE;
    }
}
