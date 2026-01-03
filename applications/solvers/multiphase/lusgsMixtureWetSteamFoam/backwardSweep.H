{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), celli)
    {
        const labelList& cellFaces = mesh.cells()[celli];

        //scalar rRho  = 0;
        //vector rRhoU = vector(0, 0, 0);
        //scalar rRhoE = 0;

        scalar rRho  = mesh.V()[celli]/dt.value()*resRho[celli];
        vector rRhoU = mesh.V()[celli]/dt.value()*resRhoU[celli];
        scalar rRhoE = mesh.V()[celli]/dt.value()*resRhoE[celli];
        scalar rRhow = mesh.V()[celli]/dt.value()*resRhow[celli];
            
        forAll(cellFaces, i)
        {
            label faceI = cellFaces[i];
            if (mesh.isInternalFace(faceI))
            {
                label nei = neighbour[faceI];
                label own = owner[faceI];

                // we use the mesh ordering property: own < nei
                if (nei != celli)  // then data[nei] belongs to upper triangle
                {
                    scalar ac = lusgsOmega*(mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar p1 = p[nei] + dp[nei];
                    vector U1 = U[nei] + dU[nei];
                    scalar T1 = T[nei] + dT[nei];
                    scalar w1 = w[nei] + dw[nei];
                    
                    scalar rho1 = 1.0/(w1/gasPropsl.rho(p1, T1) + (1.0 - w1)/gasPropsv.rho(p1, T1));
                    scalar E1 = w1*gasPropsl.Es(p1, T1) + (1.0 - w1)*gasPropsv.Es(p1, T1) + 0.5*magSqr(U1);

                    scalar phi0 = U[nei] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];
                    
                    rRho  += 0.5*ac*(rho1    - rho[nei])        - 0.5*(phi1*rho1    - phi0*rho[nei]); 
                    rRhoU += 0.5*ac*(rho1*U1 - rho[nei]*U[nei]) - 0.5*(phi1*rho1*U1 - phi0*rho[nei]*U[nei] + (p1 - p[nei])*Sf[faceI]);
                    rRhoE += 0.5*ac*(rho1*E1 - rho[nei]*E[nei]) - 0.5*(phi1*rho1*E1 - phi0*rho[nei]*E[nei] + p1*phi1 - p[nei]*phi0);
                    rRhow += 0.5*ac*(rho1*w1 - rho[nei]*w[nei]) - 0.5*(phi1*rho1*w1 - phi0*rho[nei]*w[nei]); 

                    /*scalar ac = lusgsOmega*(mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar un = U[nei] & Sf[faceI];
                    scalar dun = dU[nei] & Sf[faceI];
                    scalar unMac = un - ac;
                    
                    rRho  -= 0.5*(unMac*drhodp[nei]*dp[nei]                                                           + rho[nei]*dun                                            + unMac*drhodT[nei]*dT[nei]);
                    rRhoU -= 0.5*(unMac*drhodp[nei]*dp[nei]*U[nei]               + Sf[faceI]*dp[nei]                  + rho[nei]*dun*U[nei] + unMac*rho[nei]*dU[nei]            + unMac*drhodT[nei]*dT[nei]*U[nei]);
                    rRhoE -= 0.5*(unMac*(drhodp[nei]*H[nei] + (T[nei]/rho[nei])*drhodT[nei])*dp[nei] + un*dp[nei]     + rho[nei]*dun*H[nei] + unMac*rho[nei]*(U[nei] & dU[nei]) + unMac*(drhodT[nei]*H[nei] + rho[nei]*cp[nei])*dT[nei]);*/
                }

                if (own != celli)
                {
                    scalar ac = lusgsOmega*(mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);
                    
                    scalar p1 = p[own] + dpStar[own];
                    vector U1 = U[own] + dUStar[own];
                    scalar T1 = T[own] + dTStar[own];
                    scalar w1 = w[own] + dwStar[own];
                    
                    scalar rho1 = 1.0/(w1/gasPropsl.rho(p1, T1) + (1.0 - w1)/gasPropsv.rho(p1, T1));
                    scalar E1 = w1*gasPropsl.Es(p1, T1) + (1.0 - w1)*gasPropsv.Es(p1, T1) + 0.5*magSqr(U1);
                    
                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];
                    
                    rRho  += 0.5*ac*(rho1    - rho[own])        + 0.5*(phi1*rho1    - phi0*rho[own]); 
                    rRhoU += 0.5*ac*(rho1*U1 - rho[own]*U[own]) + 0.5*(phi1*rho1*U1 - phi0*rho[own]*U[own] + (p1 - p[own])*Sf[faceI]);
                    rRhoE += 0.5*ac*(rho1*E1 - rho[own]*E[own]) + 0.5*(phi1*rho1*E1 - phi0*rho[own]*E[own] + p1*phi1 - p[own]*phi0);
                    rRhow += 0.5*ac*(rho1*w1 - rho[own]*w[own]) + 0.5*(phi1*rho1*w1 - phi0*rho[own]*w[own]); 
                    
                    /*scalar ac = lusgsOmega*(mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);

                    scalar un = U[own] & Sf[faceI];
                    scalar dun = dU[own] & Sf[faceI];
                    scalar unMac = un + ac;
                    
                    rRho  += 0.5*(unMac*drhodp[own]*dp[own]                                                           + rho[own]*dun                                            + unMac*drhodT[own]*dT[own]);
                    rRhoU += 0.5*(unMac*drhodp[own]*dp[own]*U[own]               + Sf[faceI]*dp[own]                  + rho[own]*dun*U[own] + unMac*rho[own]*dU[own]            + unMac*drhodT[own]*dT[own]*U[own]);
                    rRhoE += 0.5*(unMac*(drhodp[own]*H[own] + (T[own]/rho[own])*drhodT[own])*dp[own] + un*dp[own]     + rho[own]*dun*H[own] + unMac*rho[own]*(U[own] & dU[own]) + unMac*(drhodT[own]*H[own] + rho[own]*cp[own])*dT[own]);*/
                }
            }
        }

        //dp[celli] = dpStar[celli] + ((subC[celli]/subA[celli])*rRho + (((drhodT[celli]*rho[celli])/subA[celli])*U[celli] & rRhoU)  + ((-drhodT[celli]*rho[celli])/subA[celli])*rRhoE)/D[celli];
        //dU[celli] = dUStar[celli] + ((-U[celli]/rho[celli])   *rRho + (1/rho[celli])*rRhoU)/D[celli];
        //dT[celli] = dTStar[celli] + ((subB[celli]/subA[celli])*rRho + (((-drhodp[celli]*rho[celli])/subA[celli])*U[celli] & rRhoU) + ((drhodp[celli]*rho[celli])/subA[celli])*rRhoE)/D[celli];
        
        //dw[celli] = dwStar[celli] + ((-w[celli]/rho[celli])   *rRho + (1.0/rho[celli])*rRhow)/D[celli];

        dp[celli] = ((invA[celli]/invD[celli])*rRho + ((drhodT[celli]/invD[celli])*U[celli] & rRhoU) + (-drhodT[celli]/invD[celli])*rRhoE + (invC[celli]/invD[celli])*rRhow)/D[celli];
        dU[celli] = ((-U[celli]/rho[celli])   *rRho + (1.0/rho[celli])*rRhoU)/D[celli];
        dT[celli] = ((invB[celli]/invD[celli])*rRho + ((-drhodp[celli]/invD[celli])*U[celli] & rRhoU) + (drhodp[celli]/invD[celli])*rRhoE + (invE[celli]/invD[celli])*rRhow)/D[celli];
        dw[celli] = ((-w[celli]/rho[celli])   *rRho + (1.0/rho[celli])*rRhow)/D[celli];
    }
}
