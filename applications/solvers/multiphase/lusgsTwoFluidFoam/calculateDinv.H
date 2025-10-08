{
    volScalarField::Internal D = mesh.V()/dt;

    forAll(owner, faceI)
    {
        const scalar ac = 0.5*lusgsOmega*( 
                            max(max(mag(U1[owner[faceI]] & Sf[faceI]) + a1[owner[faceI]]*magSf[faceI],
                                mag(U2[owner[faceI]] & Sf[faceI])     + a2[owner[faceI]]*magSf[faceI]),
                            max(mag(U1[neighbour[faceI]] & Sf[faceI]) + a1[neighbour[faceI]]*magSf[faceI],
                                mag(U2[neighbour[faceI]] & Sf[faceI]) + a2[neighbour[faceI]]*magSf[faceI])));

        D[owner[faceI]]     += ac;
        D[neighbour[faceI]] += ac;
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], faceI)
        {
            D[pFaceCells[faceI]] += 0.5*lusgsOmega*(max(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + a1[pFaceCells[faceI]]*magSf[faceI],
                                                        mag(U2[pFaceCells[faceI]] & Sf[faceI]) + a2[pFaceCells[faceI]]*magSf[faceI]));
        }
    }

    forAll(mesh.cells(), celli)
    {
        Eigen::Matrix<double, 10, 10> J = Eigen::Matrix<double, 10, 10>::Zero();

        J(0, 0) += alpha1*drhodp1[celli];
        J(1, 0) += alpha1*drhodp1[celli]*U1[celli].x();
        J(2, 0) += alpha1*drhodp1[celli]*U1[celli].y();
        J(3, 0) += alpha1*drhodp1[celli]*U1[celli].z();
        J(4, 0) += alpha1*drhodp1[celli]*H1[celli] + alpha1*drhodT1[celli]*T1[celli]/rho1[celli];
        J(5, 0) += alpha2*drhodp2[celli];
        J(6, 0) += alpha2*drhodp2[celli]*U2[celli].x();
        J(7, 0) += alpha2*drhodp2[celli]*U2[celli].y();
        J(8, 0) += alpha2*drhodp2[celli]*U2[celli].z();
        J(9, 0) += alpha2*drhodp2[celli]*H2[celli] + alpha2*drhodT2[celli]*T2[celli]/rho2[celli];

        J(0, 1) += -rho1[celli];
        J(1, 1) += -rho1[celli]*U1[celli].x();
        J(2, 1) += -rho1[celli]*U1[celli].y();
        J(3, 1) += -rho1[celli]*U1[celli].z();
        J(4, 1) += -rho1[celli]*(e1[celli] + 0.5*magSqr(U1[celli])) - pInt[celli];
        J(5, 1) +=  rho2[celli]
        J(6, 1) +=  rho2[celli]*U2[celli].x();
        J(7, 1) +=  rho2[celli]*U2[celli].y();
        J(8, 1) +=  rho2[celli]*U2[celli].z();
        J(9, 1) +=  rho2[celli]*(e2[celli] + 0.5*magSqr(U2[celli])) + pInt[celli];

        J(1, 2) += alpha1*rho1[celli];
        J(4, 2) += alpha1*rho1[celli]*U1[celli].x();
        J(2, 3) += alpha1*rho1[celli];
        J(4, 3) += alpha1*rho1[celli]*U1[celli].y();
        J(3, 4) += alpha1*rho1[celli];
        J(4, 4) += alpha1*rho1[celli]*U1[celli].z();

        J(6, 5) += alpha2*rho2[celli];
        J(9, 5) += alpha2*rho2[celli]*U2[celli].x();
        J(7, 6) += alpha2*rho2[celli];
        J(9, 6) += alpha2*rho2[celli]*U2[celli].y();
        J(8, 7) += alpha2*rho2[celli];
        J(9, 7) += alpha2*rho2[celli]*U2[celli].z();

        J(0, 8) += alpha1*drhodT1[celli];
        J(1, 8) += alpha1*drhodT1[celli]*U1[celli].x();
        J(2, 8) += alpha1*drhodT1[celli]*U1[celli].y();
        J(3, 8) += alpha1*drhodT1[celli]*U1[celli].z();
        J(4, 8) += alpha1*drhodT1[celli]*H1[celli] + alpha1*rho1[celli]*gasProps1.Cp(p[celli], T1[celli]);

        J(5, 9) += alpha2*drhodT2[celli];
        J(6, 9) += alpha2*drhodT2[celli]*U2[celli].x();
        J(7, 9) += alpha2*drhodT2[celli]*U2[celli].y();
        J(8, 9) += alpha2*drhodT2[celli]*U2[celli].z();
        J(9, 9) += alpha2*drhodT2[celli]*H2[celli] + alpha2*rho2[celli]*gasProps2.Cp(p[celli], T2[celli]);
        
        Eigen::Matrix<double, 10, 10> A = (((D[celli]*Eigen::Matrix<double, 10, 10>::Identity())) * J).inverse();
        
        Dinv_p_rho[celli]  = A(0, 0);
        Dinv_p_rhoU[celli] = vector(A(0, 1), A(0, 2), A(0, 3));
        Dinv_p_rhoE[celli] = A(0, 4);

        Dinv_U_rho[celli]  = vector(A(1, 0), A(2, 0), A(3, 0));
        Dinv_U_rhoU[celli] = tensor(A(1, 1), A(1, 2), A(1, 3), A(2, 1), A(2, 2), A(2, 3), A(3, 1), A(3, 2), A(3, 3));
        Dinv_U_rhoE[celli] = vector(A(1, 4), A(2, 4), A(3, 4));
        
        Dinv_T_rho[celli]  = A(4, 0);
        Dinv_T_rhoU[celli] = vector(A(4, 1), A(4, 2), A(4, 3));
        Dinv_T_rhoE[celli] = A(4, 4);       


        Dinv_p_alphaRho1  = A(0, 0);
        Dinv_p_alphaRhoU1 = vector(A(0, 1), A(0, 2), A(0, 3));
        Dinv_p_alphaRhoE1 = A(0, 4);
        Dinv_p_alphaRho2  = A(0, 5);
        Dinv_p_alphaRhoU2 = vector(A(0, 6), A(0, 7), A(0, 8));
        Dinv_p_alphaRhoE2 = A(0, 9);

        Dinv_alpha_alphaRho1  = A(1, 0);
        Dinv_alpha_alphaRhoU1 = vector(A(1, 1), A(1, 2), A(1, 3));
        Dinv_alpha_alphaRhoE1 = A(1, 4);
        Dinv_alpha_alphaRho2  = A(1, 5);
        Dinv_alpha_alphaRhoU2 = vector(A(1, 6), A(1, 7), A(1, 8));
        Dinv_alpha_alphaRhoE2 = A(1, 9);

        Dinv_U1_alphaRho1  = vector(A(2, 0), A(3, 0), A(4, 0));
        Dinv_U1_alphaRhoU1 = tensor(A(2, 1), A(2, 2), A(2, 3), A(3, 1), A(3, 2), A(3, 3), A(4, 1), A(4, 2), A(4, 3));
        Dinv_U1_alphaRhoE1 = vector(A(2, 4), A(3, 4), A(4, 4));
        Dinv_U1_alphaRho2  = vector(A(2, 5), A(3, 5), A(4, 5));
        Dinv_U1_alphaRhoU2 = tensor(A(2, 6), A(2, 7), A(2, 7), A(3, 6), A(3, 7), A(3, 8), A(4, 6), A(4, 7), A(4, 8));
        Dinv_U1_alphaRhoE2 = vector(A(2, 9), A(3, 9), A(4, 9));

        Dinv_U2_alphaRho1  = vector(A(5, 0), A(6, 0), A(7, 0));
        Dinv_U2_alphaRhoU1 = tensor(A(5, 1), A(5, 2), A(5, 3), A(6, 1), A(6, 2), A(6, 3), A(7, 1), A(7, 2), A(7, 3));
        Dinv_U2_alphaRhoE1 = vector(A(5, 4), A(6, 4), A(7, 4));
        Dinv_U2_alphaRho2  = vector(A(5, 5), A(6, 5), A(7, 5));
        Dinv_U2_alphaRhoU2 = tensor(A(5, 6), A(5, 7), A(5, 7), A(6, 6), A(6, 7), A(6, 8), A(7, 6), A(7, 7), A(7, 8));
        Dinv_U2_alphaRhoE2 = vector(A(5, 9), A(6, 9), A(7, 9));
        
        Dinv_T1_alphaRho1  = A(8, 0);
        Dinv_T1_alphaRhoU1 = vector(A(8, 1), A(8, 2), A(8, 3));
        Dinv_T1_alphaRhoE1 = A(8, 4);
        Dinv_T1_alphaRho2  = A(8, 5);
        Dinv_T1_alphaRhoU2 = vector(A(8, 6), A(8, 7), A(8, 8));
        Dinv_T1_alphaRhoE2 = A(8, 9);

        Dinv_T2_alphaRho1  = A(9, 0);
        Dinv_T2_alphaRhoU1 = vector(A(9, 1), A(9, 2), A(9, 3));
        Dinv_T2_alphaRhoE1 = A(9, 4);
        Dinv_T2_alphaRho2  = A(9, 5);
        Dinv_T2_alphaRhoU2 = vector(A(9, 6), A(9, 7), A(9, 8));
        Dinv_T2_alphaRhoE2 = A(9, 9);
    }
}
