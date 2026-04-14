{
    volScalarField::Internal dVdt = mesh.V()/dt;

    volScalarField::Internal acSum = 0.0*dVdt;

    //volScalarField::Internal acSum1 = 0.0*dVdt;
    //volScalarField::Internal acSum2 = 0.0*dVdt;

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

        /*acSum[owner[facei]]     += 0.5*lusgsOmega*max(max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],     mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]),
                                                      max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei], mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]));
        acSum[neighbour[facei]] += 0.5*lusgsOmega*max(max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],     mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]),
                                                      max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei], mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]));*/

        acSum[owner[facei]]     += 0.5*lusgsOmega*max(mag(U1[owner[facei]] & Sf[facei])     + a1[owner[facei]]*magSf[facei],     mag(U2[owner[facei]] & Sf[facei])     + a2[owner[facei]]*magSf[facei]);
        acSum[neighbour[facei]] += 0.5*lusgsOmega*max(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei], mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]);

        //acSum1[owner[facei]]     += 0.5*lusgsOmega*(mag(U1[owner[facei]]     & Sf[facei]) + a1[owner[facei]]    *magSf[facei]);
        //acSum2[owner[facei]]     += 0.5*lusgsOmega*(mag(U2[owner[facei]]     & Sf[facei]) + a2[owner[facei]]    *magSf[facei]);
        //acSum1[neighbour[facei]] += 0.5*lusgsOmega*(mag(U1[neighbour[facei]] & Sf[facei]) + a1[neighbour[facei]]*magSf[facei]);
        //acSum2[neighbour[facei]] += 0.5*lusgsOmega*(mag(U2[neighbour[facei]] & Sf[facei]) + a2[neighbour[facei]]*magSf[facei]);
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& pMagSf = mesh.magSf().boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            acSum[pFaceCells[facei]] += 0.5*lusgsOmega*(max(mag(U1[pFaceCells[facei]] & pSf[facei]) + a1[pFaceCells[facei]]*pMagSf[facei], mag(U2[pFaceCells[facei]] & pSf[facei]) + a2[pFaceCells[facei]]*pMagSf[facei]));

            //acSum1[pFaceCells[facei]] += 0.5*lusgsOmega*(mag(U1[pFaceCells[facei]] & pSf[facei]) + a1[pFaceCells[facei]]*pMagSf[facei]);
            //acSum2[pFaceCells[facei]] += 0.5*lusgsOmega*(mag(U2[pFaceCells[facei]] & pSf[facei]) + a2[pFaceCells[facei]]*pMagSf[facei]);
        }
    }

    forAll(mesh.cells(), celli)
    {
        //DRAG
        std::array<scalar, 100> dragJ_ = drag.dSdpUT(celli);
        Eigen::Map<Eigen::Matrix<double, 10, 10>> dragJ(dragJ_.data());

        Eigen::Matrix<double, 10, 10> virtualMassAdditionJ = Eigen::Matrix<double, 10, 10>::Zero();

        virtualMassAdditionJ(1, 2) = virtualMassCoeffs[celli];
        virtualMassAdditionJ(4, 2) = virtualMassCoeffs[celli]*U1[celli].x();
        virtualMassAdditionJ(2, 3) = virtualMassCoeffs[celli];
        virtualMassAdditionJ(4, 3) = virtualMassCoeffs[celli]*U1[celli].y();
        virtualMassAdditionJ(3, 4) = virtualMassCoeffs[celli];
        virtualMassAdditionJ(4, 4) = virtualMassCoeffs[celli]*U1[celli].z();

        virtualMassAdditionJ(6, 5) = virtualMassCoeffs[celli];
        virtualMassAdditionJ(9, 5) = virtualMassCoeffs[celli]*U2[celli].x();
        virtualMassAdditionJ(7, 6) = virtualMassCoeffs[celli];
        virtualMassAdditionJ(9, 6) = virtualMassCoeffs[celli]*U2[celli].y();
        virtualMassAdditionJ(8, 7) = virtualMassCoeffs[celli];
        virtualMassAdditionJ(9, 7) = virtualMassCoeffs[celli]*U2[celli].z();

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



        /////////////////////
        /*const labelList& cellFaces = mesh.cells()[celli];

        forAll(cellFaces,i)
        {
            label facei = cellFaces[i];

            if (!mesh.isInternalFace(facei))
            {
                label patchID = mesh.boundaryMesh().whichPatch(facei);
                label localFaceID = facei - mesh.boundaryMesh()[patchID].start();

                if(isA<wallPolyPatch>(mesh.boundaryMesh()[patchID]))
                {
                    const scalar ac1 = lusgsOmega*(max(mag(Uw1 & Sf[facei]) + (a1[celli] * magSf[facei]), mag(Uw2 & Sf[facei]) + (a2[celli]*magSf[facei])));
                    const scalar ac2 = ac1;

                    const scalar un1 = U1[celli] & Sf[facei];
                    const scalar un2 = U2[celli] & Sf[facei];
                    const scalar unMLambda1 = ac1 - un1;
                    const scalar unMLambda2 = ac2 - un2;

                    0.5*(                                                     alpha1[celli]*unMLambda1*drhodp1[celli]*dpStar[celli]                            - unMLambda1*rho1[celli]*dalphaStar[celli]                                + ((alpha1[celli]*rho1[celli]*Sf[facei]) & dUw1Star)                                + alpha1[celli]*(unMLambda1*drhodT1[celli])*dT1Star[celli]);
                    0.5*(                                                     alpha2[celli]*unMLambda2*drhodp2[celli]*dpStar[celli]                            + unMLambda2*rho2[celli]*dalphaStar[celli]                                + ((alpha2[celli]*rho2[celli]*Sf[facei]) & dUw2Star)                                + alpha2[celli]*(unMLambda2*drhodT2[celli])*dT2Star[celli]);

                    0.5*(                                   alpha1[celli]*(unMLambda1*drhodp1[celli]*Uw1 + Sf[facei])*dpStar[celli] - (unMLambda1*rho1[celli]*Uw1 + p[celli]*Sf[facei])*dalphaStar[celli] + (alpha1[celli]*rho1[celli]*(unMLambda1*tensor::I + (Uw1 * Sf[facei])) & dUw1Star)                            + alpha1[celli]*(unMLambda1*drhodT1[celli]*Uw1)*dT1Star[celli]);
                    0.5*(                                   alpha2[celli]*(unMLambda2*drhodp2[celli]*Uw2 + Sf[facei])*dpStar[celli] + (unMLambda2*rho2[celli]*Uw2 + p[celli]*Sf[facei])*dalphaStar[celli] + (alpha2[celli]*rho2[celli]*(unMLambda2*tensor::I + (Uw2 * Sf[facei])) & dUw2Star)                            + alpha2[celli]*(unMLambda2*drhodT2[celli]*Uw2)*dT2Star[celli]);

                    0.5*(alpha1[celli]*(unMLambda1*(drhodp1[celli]*Hw1 + drhodT1[celli]*T1[celli]/rho1[celli]) + un1)*dpStar[celli]       - (unMLambda1*rho1[celli]*Ew1 + p[celli]*un1)*dalphaStar[celli]         + ((alpha1[celli]*rho1[celli]*(unMLambda1*Uw1 + Hw1*Sf[facei])) & dUw1Star) + alpha1[celli]*(unMLambda1*(drhodT1[celli]*Hw1 + rho1[celli]*cp1[celli]))*dT1Star[celli]);
                    0.5*(alpha2[celli]*(unMLambda2*(drhodp2[celli]*Hw2 + drhodT2[celli]*T2[celli]/rho2[celli]) + un2)*dpStar[celli]       + (unMLambda2*rho2[celli]*Ew2 + p[celli]*un2)*dalphaStar[celli]         + ((alpha2[celli]*rho2[celli]*(unMLambda2*Uw2 + Hw2*Sf[facei])) & dUw2Star) + alpha2[celli]*(unMLambda2*(drhodT2[celli]*Hw2 + rho2[celli]*cp2[celli]))*dT2Star[celli]);


                    J(0, 0) += 0.5*alpha1[celli]*(un1*drhodp1[celli]);
                    J(1, 0) += 0.5*alpha1[celli]*(un1*drhodp1[celli]*U1.x() + Sf[facei].x())
                    J(2, 0) += 0.5*alpha1[celli]*(un1*drhodp1[celli]*U1.y() + Sf[facei].y())
                    J(3, 0) += 0.5*alpha1[celli]*(un1*drhodp1[celli]*U1.z() + Sf[facei].z())
                    J(4, 0) += 0.5*alpha1[celli]*(un1*(drhodp1[celli]*H1[celli] + drhodT1[celli]*T1[celli]/rho1[celli]) + un1);
                    J(5, 0) += 0.5*alpha2[celli]*(un2*drhodp2[celli]);
                    J(6, 0) += 0.5*alpha2[celli]*(un2*drhodp2[celli]*U2.x() + Sf[facei].x())
                    J(7, 0) += 0.5*alpha2[celli]*(un2*drhodp2[celli]*U2.y() + Sf[facei].y())
                    J(8, 0) += 0.5*alpha2[celli]*(un2*drhodp2[celli]*U2.z() + Sf[facei].z())
                    J(9, 0) += 0.5*alpha2[celli]*(un2*(drhodp2[celli]*H2[celli] + drhodT2[celli]*T2[celli]/rho2[celli]) + un2);

                    J(0, 1) += -0.5*un1*rho1[celli];
                    J(1, 1) += -0.5*(un1*rho1[celli]*U1[celli].x() + p[celli]*Sf[facei].x());
                    J(2, 1) += -0.5*(un1*rho1[celli]*U1[celli].y() + p[celli]*Sf[facei].y());
                    J(3, 1) += -0.5*(un1*rho1[celli]*U1[celli].z() + p[celli]*Sf[facei].z());
                    J(4, 1) += -0.5*(un1*rho1[celli]*(e1[celli] + 0.5*magSqr(U1[celli])) + p[celli]*un1);
                    J(5, 1) +=  0.5*un2*rho2[celli];
                    J(6, 1) +=  0.5*(un2*rho2[celli]*U2[celli].x() + p[celli]*Sf[facei].x());
                    J(7, 1) +=  0.5*(un2*rho2[celli]*U2[celli].y() + p[celli]*Sf[facei].y());
                    J(8, 1) +=  0.5*(un2*rho2[celli]*U2[celli].z() + p[celli]*Sf[facei].z());
                    J(9, 1) +=  0.5*(un2*rho2[celli]*(e2[celli] + 0.5*magSqr(U2[celli])) + p[celli]*un2);


                    //TODO ->

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
                }
            }

        }*/
        ////////////////



        D_inv[celli].setThreshold(1e-60);     
        //Eigen::Matrix<double, 10, 10> A = (D[celli]*Eigen::Matrix<double, 10, 10>::Identity()) * J;
        //D_inv[celli].compute(A);
   
        //D_inv[celli].compute(dVdt[celli] * (J + pIntAdditionJ + virtualMassAdditionJ) + ((Eigen::Vector<double, 10>() << Eigen::Vector<double, 5>::Constant(acSum1[celli]), Eigen::Vector<double, 5>::Constant(acSum2[celli])).finished().asDiagonal() * J) + (dragJ*mesh.V()[celli]));
        D_inv[celli].compute(dVdt[celli] * (J + pIntAdditionJ + virtualMassAdditionJ) + (Eigen::Vector<double, 10>::Constant(acSum[celli]).asDiagonal() * J) + (dragJ*mesh.V()[celli]));

        //alpha1zero[celli] = (alpha1[celli] < 1.0e-20) ? 0.0 : 1.0;
        //alpha2zero[celli] = (alpha2[celli] < 1.0e-20) ? 0.0 : 1.0;

        //SVD
        //D_inv[celli]->setThreshold(1e-20);
        //D_inv[celli]->compute(dVdt[celli] * (J + pIntAdditionJ) + (Eigen::Vector<double, 10>::Constant(acSum[celli]).asDiagonal() * J) + (dragJ*mesh.V()[celli]), Eigen::ComputeFullU | Eigen::ComputeFullV);
    }
}
