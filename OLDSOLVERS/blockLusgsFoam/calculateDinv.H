{
    volScalarField::Internal D = mesh.V()/dt;

    forAll(owner, faceI)
    {
        D[owner[faceI]]     += 0.5*lusgsOmega*(mag(U[owner[faceI]]     & Sf[faceI]) + a[owner[faceI]]    *magSf[faceI]);
        D[neighbour[faceI]] += 0.5*lusgsOmega*(mag(U[neighbour[faceI]] & Sf[faceI]) + a[neighbour[faceI]]*magSf[faceI]);
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], faceI)
        {
            D[pFaceCells[faceI]] += 0.5*lusgsOmega*(mag(U[pFaceCells[faceI]] & Sf[faceI]) + a[pFaceCells[faceI]]*magSf[faceI]);
        }
    }

    forAll(mesh.cells(), celli)
    {
        Eigen::Matrix<double, 5, 5> J = Eigen::Matrix<double, 5, 5>::Zero();

        J(0, 0) += drhodp[celli];
        J(1, 0) += drhodp[celli]*U[celli].x();
        J(2, 0) += drhodp[celli]*U[celli].y();
        J(3, 0) += drhodp[celli]*U[celli].z();
        J(4, 0) += drhodp[celli]*H[celli] + drhodT[celli]*T[celli]/rho[celli];

        J(1, 1) += rho[celli];
        J(4, 1) += rho[celli]*U[celli].x();
        J(2, 2) += rho[celli];
        J(4, 2) += rho[celli]*U[celli].y();
        J(3, 3) += rho[celli];
        J(4, 3) += rho[celli]*U[celli].z();

        J(0, 4) += drhodT[celli];
        J(1, 4) += drhodT[celli]*U[celli].x();
        J(2, 4) += drhodT[celli]*U[celli].y();
        J(3, 4) += drhodT[celli]*U[celli].z();
        J(4, 4) += drhodT[celli]*H[celli] + rho[celli]*gasProps.Cp(p[celli], T[celli]);

        //Eigen::Matrix<double, 5, 5> A = ((D[celli]*Eigen::Matrix<double, 5, 5>::Identity())*J);
        
        Eigen::Matrix<double, 5, 5> A = (((D[celli]*Eigen::Matrix<double, 5, 5>::Identity())) * J).inverse();
        
        Dinv_p_rho[celli]  = A(0, 0);
        Dinv_p_rhoU[celli] = vector(A(0, 1), A(0, 2), A(0, 3));
        Dinv_p_rhoE[celli] = A(0, 4);

        Dinv_U_rho[celli]  = vector(A(1, 0), A(2, 0), A(3, 0));
        Dinv_U_rhoU[celli] = tensor(A(1, 1), A(1, 2), A(1, 3), A(2, 1), A(2, 2), A(2, 3), A(3, 1), A(3, 2), A(3, 3));
        Dinv_U_rhoE[celli] = vector(A(1, 4), A(2, 4), A(3, 4));
        
        Dinv_T_rho[celli]  = A(4, 0);
        Dinv_T_rhoU[celli] = vector(A(4, 1), A(4, 2), A(4, 3));
        Dinv_T_rhoE[celli] = A(4, 4);


        /*pCond_rho_p[celli]  = J(0, 0);
        pCond_rho_U[celli] = vector(J(0, 1), J(0, 2), J(0, 3));
        pCond_rho_T[celli] = J(0, 4);

        pCond_rhoU_p[celli]  = vector(J(1, 0), J(2, 0), J(3, 0));
        pCond_rhoU_U[celli] = tensor(J(1, 1), J(1, 2), J(1, 3), J(2, 1), J(2, 2), J(2, 3), J(3, 1), J(3, 2), J(3, 3));
        pCond_rhoU_T[celli] = vector(J(1, 4), J(2, 4), J(3, 4));
        
        pCond_rhoE_p[celli]  = J(4, 0);
        pCond_rhoE_U[celli] = vector(J(4, 1), J(4, 2), J(4, 3));
        pCond_rhoE_T[celli] = J(4, 4);*/
        
    }
}
