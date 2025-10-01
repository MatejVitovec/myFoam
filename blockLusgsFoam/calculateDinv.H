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
        const scalar rho_ = thermo.rho()[celli];

        const scalar drhodp = gasProps.beta_T(p[celli], T[celli]);
        const scalar drhodT = gasProps.beta_p(p[celli], T[celli]);
        const scalar cp = gasProps.Cp(p[celli], T[celli]);

        const scalar H = thermo.he()[celli] + 0.5*magSqr(U[celli]) + p[celli]/rho;

        Eigen::Matrix<double, 5, 5> J = Eigen::Matrix<double, 5, 5>::Zero();

        J(0, 0) += drhodp;
        J(1, 0) += drhodp*U[celli].x();
        J(2, 0) += drhodp*U[celli].y();
        J(3, 0) += drhodp*U[celli].z();
        J(4, 0) += drhodp*H + drhodT*T[celli]/rho_;

        J(1, 1) += U[celli].x();
        J(4, 1) += U[celli].x()*rho_;
        J(2, 2) += U[celli].y();
        J(4, 2) += U[celli].y()*rho_;
        J(3, 3) += U[celli].z();
        J(4, 3) += U[celli].z()*rho_;

        J(0, 4) += drhodT;
        J(1, 4) += drhodT*U[celli].x();;
        J(2, 4) += drhodT*U[celli].z();;
        J(3, 4) += drhodT*U[celli].y();;
        J(4, 4) += drhodT*H + rho_*gasProps.Cp(p[celli], T[celli]);

        //Eigen::Matrix<double, 5, 5> A = ((D[celli]*Eigen::Matrix<double, 5, 5>::Identity())*J);
        
        Eigen::Matrix<double, 5, 5> A = ((D[celli]*Eigen::Matrix<double, 5, 5>::Identity())*J).inverse();

        
        Dinv_p_rho[celli]  = A(0, 0);
        Dinv_p_rhoU[celli] = vector(A(0, 1), A(0, 2), A(0, 3));
        Dinv_p_rhoE[celli] = A(0, 4);

        Dinv_U_rho[celli]  = vector(A(1, 0), A(2, 0), A(3, 0));
        Dinv_U_rhoU[celli] = tensor(A(1, 1), A(1, 2), A(1, 3), A(2, 1), A(2, 2), A(2, 3), A(3, 1), A(3, 2), A(3, 3));
        Dinv_U_rhoE[celli] = vector(A(1, 4), A(2, 4), A(3, 4));
        
        Dinv_T_rho[celli]  = A(4, 0);
        Dinv_T_rhoU[celli] = vector(A(4, 1), A(4, 2), A(4, 3));
        Dinv_T_rhoE[celli] = A(4, 4);
    }

}