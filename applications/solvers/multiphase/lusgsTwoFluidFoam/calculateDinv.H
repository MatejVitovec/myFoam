{
    volScalarField::Internal dVdt = mesh.V()/dt;

    volScalarField::Internal acSum = 0.0*dVdt;

    volScalarField::Internal D1 = 0.0*dVdt;
    volScalarField::Internal D2 = 0.0*dVdt;

    volScalarField::Internal ac1 = 0.0*dVdt;
    volScalarField::Internal ac2 = 0.0*dVdt;
    
    forAll(owner, facei)
    {
        /*const scalar ac = 0.5*lusgsOmega*( 
                            max(max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],
                                    mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]),
                                max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei],
                                    mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei])));*/

        //acSum[owner[facei]]     += 0.5*lusgsOmega*max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],     mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]);
        //acSum[neighbour[facei]] += 0.5*lusgsOmega*max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei], mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]);

        acSum[owner[facei]]     += 0.5*ac[facei];
        acSum[neighbour[facei]] += 0.5*ac[facei];

        //D1[owner[facei]]     += 0.5*lusgsOmega*(mag(U1[owner[facei]]     & Sf[facei]) + a1[owner[facei]]    *magSf[facei]);
        //D2[owner[facei]]     += 0.5*lusgsOmega*(mag(U2[owner[facei]]     & Sf[facei]) + a2[owner[facei]]    *magSf[facei]);
        //D1[neighbour[facei]] += 0.5*lusgsOmega*(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei]);
        //D2[neighbour[facei]] += 0.5*lusgsOmega*(mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], facei)
        {
            //acSum[pFaceCells[facei]] += 0.5*lusgsOmega*(max(mag(U1[pFaceCells[facei]] & Sf[facei]) + a1[pFaceCells[facei]]*magSf[facei], mag(U2[pFaceCells[facei]] & Sf[facei]) + a2[pFaceCells[facei]]*magSf[facei]));
            //acSum[pFaceCells[facei]] += 0.5*ac[facei];

            //D1[pFaceCells[facei]] += 0.5*lusgsOmega*(mag(U1[pFaceCells[facei]] & Sf[facei]) + a1[pFaceCells[facei]]*magSf[facei]);
            //D2[pFaceCells[facei]] += 0.5*lusgsOmega*(mag(U2[pFaceCells[facei]] & Sf[facei]) + a2[pFaceCells[facei]]*magSf[facei]);
        }
    }

    forAll(mesh.cells(), celli)
    {
        //DRAG
        const double K = dragK[celli];

        std::array<scalar, 100> dragJ_ = drag.dSdpUT(celli, K, d[celli]);
        Eigen::Map<Eigen::Matrix<double, 10, 10>> dragJ(dragJ_.data());

        /*Eigen::Matrix<double, 10, 10> dragJ = Eigen::Matrix<double, 10, 10>::Zero();
        
        dragJ(1, 2) = K;
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
        //END DRAG

        Eigen::Matrix<double, 10, 10> J = Eigen::Matrix<double, 10, 10>::Zero();
        Eigen::Matrix<double, 10, 10> pIntAdditionJ = Eigen::Matrix<double, 10, 10>::Zero();

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
        J(4, 1) += -rho1[celli]*(e1[celli] + 0.5*magSqr(U1[celli])); //- fluid.pInt()[celli];
        J(5, 1) +=  rho2[celli];
        J(6, 1) +=  rho2[celli]*U2[celli].x();
        J(7, 1) +=  rho2[celli]*U2[celli].y();
        J(8, 1) +=  rho2[celli]*U2[celli].z();
        J(9, 1) +=  rho2[celli]*(e2[celli] + 0.5*magSqr(U2[celli])); //+ fluid.pInt()[celli];

        pIntAdditionJ(4, 1) += -fluid.pInt()[celli];
        pIntAdditionJ(9, 1) +=  fluid.pInt()[celli];

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
        //D_inv[celli].compute(((Eigen::Vector<double, 10>() << Eigen::Vector<double, 5>::Constant(D1[celli]), Eigen::Vector<double, 5>::Constant(D2[celli])).finished().asDiagonal() * J) + dragJ);
        
        //D_inv[celli].compute(dVdt[celli] * (J + pIntAdditionJ) + ((Eigen::Vector<double, 10>() << Eigen::Vector<double, 5>::Constant(D1[celli]), Eigen::Vector<double, 5>::Constant(D2[celli])).finished().asDiagonal() * J) + (dragJ*mesh.V()[celli]));
        D_inv[celli].compute(dVdt[celli] * (J + pIntAdditionJ) + (Eigen::Vector<double, 10>::Constant(acSum[celli]).asDiagonal() * J) + (dragJ*mesh.V()[celli]));

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
