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
                    scalar ac = 0.5*lusgsOmega*( 
                        mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);

                    scalar rhoE = rho[own]*(e[own] + 0.5*magSqr(U[own]));

                    scalar p1 = p[own] + dpStar[own];
                    vector U1 = U[own] + dUStar[own];
                    scalar T1 = T[own] + dTStar[own];
                    
                    scalar rho1 = gasProps.rho(p1, T1);
                    scalar rhoE1 = rho1*(gasProps.Es(p1, T1) + 0.5*magSqr(U1));

                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];

                    /*scalar deltaRho  = pCond_rho_p[cellI] *dpStar[cellI] + (pCond_rho_U[cellI]  & dUStar[cellI]) + pCond_rho_T[cellI] *dTStar[cellI];
                    vector deltaRhoU = pCond_rhoU_p[cellI]*dpStar[cellI] + (pCond_rhoU_U[cellI] & dUStar[cellI]) + pCond_rhoU_T[cellI]*dTStar[cellI];
                    scalar deltaRhoE = pCond_rhoE_p[cellI]*dpStar[cellI] + (pCond_rhoE_U[cellI] & dUStar[cellI]) + pCond_rhoE_T[cellI]*dTStar[cellI];

                    rRho  += ac*deltaRho  - 0.5*(rho1*phi1 - rho[own]*phi0); 
                    rRhoU += ac*deltaRhoU - 0.5*((rho1*U1*phi1 - rho[own]*U[own]*phi0) + (p1 - p[own])*Sf[faceI]);
                    rRhoE += ac*deltaRhoE - 0.5*(rhoE1*phi1 + p1*phi1 - (rhoE*phi0 + p[own]*phi0)); */

                    rRho  += ac*(rho1    - rho[own])        - 0.5*(rho1*phi1 - rho[own]*phi0); 
                    rRhoU += ac*(rho1*U1 - rho[own]*U[own]) - 0.5*((rho1*U1*phi1 - rho[own]*U[own]*phi0) + (p1 - p[own])*Sf[faceI]);
                    rRhoE += ac*(rhoE1   - rhoE)            - 0.5*(rhoE1*phi1 + p1*phi1 - (rhoE*phi0 + p[own]*phi0)); 
                }

                /*if (nei != cellI)  // then data[own] belongs to lower triangle
                {
                    scalar ac = 0.5*lusgsOmega*( 
                        mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar rhoE = rho[nei]*(e[nei] + 0.5*magSqr(U[nei]));

                    scalar p1 = p[nei] + dp[nei];
                    vector U1 = U[nei] + dU[nei];
                    scalar T1 = T[nei] + dT[nei];
                    
                    scalar rho1 = gasProps.rho(p1, T1);
                    scalar rhoE1 = rho1*(gasProps.Es(p1, T1) + 0.5*magSqr(U1));

                    scalar phi0 = U[nei] & Sf[faceI];
                    scalar phi1 = U1 & Sf[faceI];

                    rRho  += ac*(rho1    - rho[nei])        - 0.5*(rho1*phi1 - rho[nei]*phi0); 
                    rRhoU += ac*(rho1*U1 - rho[nei]*U[nei]) - 0.5*((rho1*U1*phi1 - rho[nei]*U[nei]*phi0) + (p1 - p[nei])*Sf[faceI]);
                    rRhoE += ac*(rhoE1   - rhoE)            - 0.5*(rhoE1*phi1 + p1*phi1 - (rhoE*phi0 + p[nei]*phi0)); 
                }*/
            }
        }

        dpStar[cellI] = Dinv_p_rho[cellI]*rRho + (Dinv_p_rhoU[cellI] & rRhoU) + Dinv_p_rhoE[cellI]*rRhoE;
        dUStar[cellI] = Dinv_U_rho[cellI]*rRho + (Dinv_U_rhoU[cellI] & rRhoU) + Dinv_U_rhoE[cellI]*rRhoE;
        dTStar[cellI] = Dinv_T_rho[cellI]*rRho + (Dinv_T_rhoU[cellI] & rRhoU) + Dinv_T_rhoE[cellI]*rRhoE;
    }
}
