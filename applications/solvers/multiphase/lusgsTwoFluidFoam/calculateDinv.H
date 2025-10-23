{
    //volScalarField::Internal D = mesh.V()/dt;

    volScalarField::Internal D1 = mesh.V()/dt;
    volScalarField::Internal D2 = mesh.V()/dt;

    forAll(owner, faceI)
    {
        /*const scalar ac = 0.5*lusgsOmega*( 
                            max(max(mag(U1[owner[faceI]] & Sf[faceI])     + a1[owner[faceI]]*magSf[faceI],
                                    mag(U2[owner[faceI]] & Sf[faceI])     + a2[owner[faceI]]*magSf[faceI]),
                                max(mag(U1[neighbour[faceI]] & Sf[faceI]) + a1[neighbour[faceI]]*magSf[faceI],
                                    mag(U2[neighbour[faceI]] & Sf[faceI]) + a2[neighbour[faceI]]*magSf[faceI])));*/

        //D[owner[faceI]]     += 0.5*lusgsOmega*max(mag(U1[owner[faceI]] & Sf[faceI])     + a1[owner[faceI]]*magSf[faceI],     mag(U2[owner[faceI]] & Sf[faceI])     + a2[owner[faceI]]*magSf[faceI]);
        //D[neighbour[faceI]] += 0.5*lusgsOmega*max(mag(U1[neighbour[faceI]] & Sf[faceI]) + a1[neighbour[faceI]]*magSf[faceI], mag(U2[neighbour[faceI]] & Sf[faceI]) + a2[neighbour[faceI]]*magSf[faceI]);

        D1[owner[faceI]]     += 0.5*lusgsOmega*(mag(U1[owner[faceI]]     & Sf[faceI]) + a1[owner[faceI]]    *magSf[faceI]);
        D2[owner[faceI]]     += 0.5*lusgsOmega*(mag(U2[owner[faceI]]     & Sf[faceI]) + a2[owner[faceI]]    *magSf[faceI]);
        D1[neighbour[faceI]] += 0.5*lusgsOmega*(mag(U1[neighbour[faceI]] & Sf[faceI]) + a1[neighbour[faceI]]*magSf[faceI]);
        D2[neighbour[faceI]] += 0.5*lusgsOmega*(mag(U2[neighbour[faceI]] & Sf[faceI]) + a2[neighbour[faceI]]*magSf[faceI]);
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], faceI)
        {
            //D[pFaceCells[faceI]] += 0.5*lusgsOmega*(max(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + a1[pFaceCells[faceI]]*magSf[faceI], mag(U2[pFaceCells[faceI]] & Sf[faceI]) + a2[pFaceCells[faceI]]*magSf[faceI]));

            D1[pFaceCells[faceI]] += 0.5*lusgsOmega*(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + a1[pFaceCells[faceI]]*magSf[faceI]);
            D2[pFaceCells[faceI]] += 0.5*lusgsOmega*(mag(U2[pFaceCells[faceI]] & Sf[faceI]) + a2[pFaceCells[faceI]]*magSf[faceI]);
        }
    }

    forAll(mesh.cells(), celli)
    {
        //DRAG
        const double K = dragK[celli];

        Eigen::Matrix<double, 10, 10> dragJ = Eigen::Matrix<double, 10, 10>::Zero();
        

        /*dragJ(1, 2) = K;
        dragJ(4, 2) = K*(2*U1[celli].x() - U2[celli].x());
        dragJ(6, 2) = -K;
        dragJ(9, 2) = -K*U2[celli].x();

        dragJ(2, 3) = K;
        dragJ(4, 3) = K*(2*U1[celli].y() - U2[celli].y());
        dragJ(7, 3) = -K;
        dragJ(9, 3) = -K*U2[celli].y();

        dragJ(3, 4) = K;
        dragJ(4, 4) = K*(2*U1[celli].z() - U2[celli].z());
        dragJ(8, 4) = -K;
        dragJ(9, 4) = -K*U2[celli].z();


        dragJ(1, 5) = -K;
        dragJ(4, 5) = -K*U1[celli].x(); 
        dragJ(6, 5) = K;
        dragJ(9, 5) = K*(2*U2[celli].x() - U1[celli].x());

        dragJ(2, 6) = -K;
        dragJ(4, 6) = -K*U1[celli].y(); 
        dragJ(7, 6) = K;
        dragJ(9, 6) = K*(2*U2[celli].y() - U1[celli].y());

        dragJ(3, 7) = -K;
        dragJ(4, 7) = -K*U1[celli].z(); 
        dragJ(8, 7) = K;
        dragJ(9, 7) = K*(2*U2[celli].z() - U1[celli].z());*/

        double A = 0;
        double D = mag(U1[celli] - U2[celli]);
        double D1 = (U1[celli] - U2[celli]) & U1[celli];
        double D2 = (U2[celli] - U1[celli]) & U2[celli];
        vector d = U1[celli] - U2[celli];

        dragJ(1, 2) = (d.x()*d.x())/D + D;
        dragJ(2, 2) = (d.x()*d.y())/D;
        dragJ(3, 2) = (d.x()*d.z())/D;
        dragJ(4, 2) = ((d.x()*D1)/D + (2*U1[celli].x() - U2[celli].x())*D);
        dragJ(6, 2) = -dragJ(1, 2);
        dragJ(7, 2) = -dragJ(2, 2);
        dragJ(8, 2) = -dragJ(3, 2);
        dragJ(9, 2) = -(-(d.x()*D2)/D + U2[celli].x()*D);

        //END DRAG

        Eigen::Matrix<double, 10, 10> J = Eigen::Matrix<double, 10, 10>::Zero();

        J(0, 0) += alpha1[celli]*drhodp1[celli];
        J(1, 0) += alpha1[celli]*drhodp1[celli]*U1[celli].x();
        J(2, 0) += alpha1[celli]*drhodp1[celli]*U1[celli].y();
        J(3, 0) += alpha1[celli]*drhodp1[celli]*U1[celli].z();
        J(4, 0) += alpha1[celli]*drhodp1[celli]*H1[celli] + alpha1[celli]*drhodT1[celli]*T1[celli]/rho1[celli];
        J(5, 0) += alpha2[celli]*drhodp2[celli];
        J(6, 0) += alpha2[celli]*drhodp2[celli]*U2[celli].x();
        J(7, 0) += alpha2[celli]*drhodp2[celli]*U2[celli].y();
        J(8, 0) += alpha2[celli]*drhodp2[celli]*U2[celli].z();
        J(9, 0) += alpha2[celli]*drhodp2[celli]*H2[celli] + alpha2[celli]*drhodT2[celli]*T2[celli]/rho2[celli];

        J(0, 1) += -rho1[celli];
        J(1, 1) += -rho1[celli]*U1[celli].x();
        J(2, 1) += -rho1[celli]*U1[celli].y();
        J(3, 1) += -rho1[celli]*U1[celli].z();
        J(4, 1) += -rho1[celli]*(e1[celli] + 0.5*magSqr(U1[celli])) - fluid.pInt()[celli];
        J(5, 1) +=  rho2[celli];
        J(6, 1) +=  rho2[celli]*U2[celli].x();
        J(7, 1) +=  rho2[celli]*U2[celli].y();
        J(8, 1) +=  rho2[celli]*U2[celli].z();
        J(9, 1) +=  rho2[celli]*(e2[celli] + 0.5*magSqr(U2[celli])) + fluid.pInt()[celli];

        J(1, 2) += alpha1[celli]*rho1[celli];
        J(4, 2) += alpha1[celli]*rho1[celli]*U1[celli].x();
        J(2, 3) += alpha1[celli]*rho1[celli];
        J(4, 3) += alpha1[celli]*rho1[celli]*U1[celli].y();
        J(3, 4) += alpha1[celli]*rho1[celli];
        J(4, 4) += alpha1[celli]*rho1[celli]*U1[celli].z();

        J(6, 5) += alpha2[celli]*rho2[celli];
        J(9, 5) += alpha2[celli]*rho2[celli]*U2[celli].x();
        J(7, 6) += alpha2[celli]*rho2[celli];
        J(9, 6) += alpha2[celli]*rho2[celli]*U2[celli].y();
        J(8, 7) += alpha2[celli]*rho2[celli];
        J(9, 7) += alpha2[celli]*rho2[celli]*U2[celli].z();

        J(0, 8) += alpha1[celli]*drhodT1[celli];
        J(1, 8) += alpha1[celli]*drhodT1[celli]*U1[celli].x();
        J(2, 8) += alpha1[celli]*drhodT1[celli]*U1[celli].y();
        J(3, 8) += alpha1[celli]*drhodT1[celli]*U1[celli].z();
        J(4, 8) += alpha1[celli]*drhodT1[celli]*H1[celli] + alpha1[celli]*rho1[celli]*gasProps1.Cp(p[celli], T1[celli]);

        J(5, 9) += alpha2[celli]*drhodT2[celli];
        J(6, 9) += alpha2[celli]*drhodT2[celli]*U2[celli].x();
        J(7, 9) += alpha2[celli]*drhodT2[celli]*U2[celli].y();
        J(8, 9) += alpha2[celli]*drhodT2[celli]*U2[celli].z();
        J(9, 9) += alpha2[celli]*drhodT2[celli]*H2[celli] + alpha2[celli]*rho2[celli]*gasProps2.Cp(p[celli], T2[celli]);

        //Eigen::Matrix<double, 10, 10> A = (D[celli]*Eigen::Matrix<double, 10, 10>::Identity()) * J;
        //D_inv[celli].compute(A);

        //D_inv[celli].compute((D[celli]*Eigen::Matrix<double, 10, 10>::Identity()) * J);

        D_inv[celli].setThreshold(1e-300);
        //Eigen::Matrix<long double, 10, 10> A = (Eigen::Vector<long double, 10>() << Eigen::Vector<long double, 5>::Constant(D1[celli]), Eigen::Vector<long double, 5>::Constant(D2[celli])).finished().asDiagonal() * J;
        D_inv[celli].compute(((Eigen::Vector<double, 10>() << Eigen::Vector<double, 5>::Constant(D1[celli]), Eigen::Vector<double, 5>::Constant(D2[celli])).finished().asDiagonal() * J) + dragJ);

        //if (runTime.timeIndex() > 8000)
        //{
            //if (celli > 498)
            //{
            /*    std::cout << A << std::endl << std::endl;

                Eigen::Matrix<long double, 10, 10> aux = D_inv[celli].inverse();
                std::cout << aux << std::endl << std::endl;

                std::cout << "alpha: " << alpha2[celli] << std::endl;
                std::cout << "D1: " << D1[celli] << std::endl;
                std::cout << "D2: " << D2[celli] << std::endl << std::endl;

                std::cin.ignore();*/
            //}
        //}
    }
}
