{
    volScalarField::Internal dVdt = mesh.V()/dt;

    //volScalarField::Internal acSum = 0.0*dVdt;

    volScalarField::Internal acSum1 = 0.0*dVdt;
    volScalarField::Internal acSum2 = 0.0*dVdt;

    volScalarField::Internal ac1 = 0.0*dVdt;
    volScalarField::Internal ac2 = 0.0*dVdt;
    
    forAll(owner, facei)
    {
        /*const scalar ac = 0.5*lusgsOmega*( 
                            max(max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],
                                    mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]),
                                max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei],
                                    mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei])));*/

        //acSum[owner[facei]]     += 0.5*lusgsOmega*(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei]);
        //acSum[neighbour[facei]] += 0.5*lusgsOmega*(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei]);
        //acSum[owner[facei]]     += 0.5*ac[facei];
        //acSum[neighbour[facei]] += 0.5*ac[facei];
        //-----------

        //acSum[owner[facei]]     += 0.5*lusgsOmega*max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],     mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]);
        //acSum[neighbour[facei]] += 0.5*lusgsOmega*max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei], mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]);

        acSum1[owner[facei]]     += 0.5*lusgsOmega*(mag(U1[owner[facei]]     & Sf[facei]) + a1[owner[facei]]    *magSf[facei]);
        acSum2[owner[facei]]     += 0.5*lusgsOmega*(mag(U2[owner[facei]]     & Sf[facei]) + a2[owner[facei]]    *magSf[facei]);
        acSum1[neighbour[facei]] += 0.5*lusgsOmega*(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei]);
        acSum2[neighbour[facei]] += 0.5*lusgsOmega*(mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            //acSum[pFaceCells[facei]] += 0.5*lusgsOmega*(max(mag(U1[pFaceCells[facei]] & pSf[facei]) + a1[pFaceCells[facei]]*pMagSf[facei], mag(U2[pFaceCells[facei]] & pSf[facei]) + a2[pFaceCells[facei]]*pMagSf[facei]));

            acSum1[pFaceCells[facei]] += 0.5*lusgsOmega*(mag(U1[pFaceCells[facei]] & pSf[facei]) + a1[pFaceCells[facei]]*pMagSf[facei]);
            acSum2[pFaceCells[facei]] += 0.5*lusgsOmega*(mag(U2[pFaceCells[facei]] & pSf[facei]) + a2[pFaceCells[facei]]*pMagSf[facei]);
        }
    }

    forAll(mesh.cells(), celli)
    {
        //DRAG
        const double Ki = dragKi[celli];
        std::array<scalar, 100> dragJ_ = drag.dSdpUT(celli, Ki, d[celli]);
        Eigen::Map<Eigen::Matrix<double, 10, 10>> dragJ(dragJ_.data());

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
        J(4, 1) += -rho1[celli]*(e1[celli] + 0.5*magSqr(U1[celli]));
        J(5, 1) +=  rho2[celli];
        J(6, 1) +=  rho2[celli]*U2[celli].x();
        J(7, 1) +=  rho2[celli]*U2[celli].y();
        J(8, 1) +=  rho2[celli]*U2[celli].z();
        J(9, 1) +=  rho2[celli]*(e2[celli] + 0.5*magSqr(U2[celli]));

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


        D_inv[celli].setThreshold(1e-60);     
        //Eigen::Matrix<double, 10, 10> A = (D[celli]*Eigen::Matrix<double, 10, 10>::Identity()) * J;
        //D_inv[celli].compute(A);
   
        D_inv[celli].compute(dVdt[celli] * (J + pIntAdditionJ) + ((Eigen::Vector<double, 10>() << Eigen::Vector<double, 5>::Constant(acSum1[celli]), Eigen::Vector<double, 5>::Constant(acSum2[celli])).finished().asDiagonal() * J) + (dragJ*mesh.V()[celli]));
        //D_inv[celli].compute(dVdt[celli] * (J + pIntAdditionJ) + (Eigen::Vector<double, 10>::Constant(acSum[celli]).asDiagonal() * J) + (dragJ*mesh.V()[celli]));

        //alpha1zero[celli] = (alpha1[celli] < 1.0e-20) ? 0.0 : 1.0;
        //alpha2zero[celli] = (alpha2[celli] < 1.0e-20) ? 0.0 : 1.0;

        //SVD
        //D_inv[celli]->setThreshold(1e-20);
        //D_inv[celli]->compute(dVdt[celli] * (J + pIntAdditionJ) + (Eigen::Vector<double, 10>::Constant(acSum[celli]).asDiagonal() * J) + (dragJ*mesh.V()[celli]), Eigen::ComputeFullU | Eigen::ComputeFullV);
    }
}
