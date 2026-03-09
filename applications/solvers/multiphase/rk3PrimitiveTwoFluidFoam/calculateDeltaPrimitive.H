{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    const gasProperties& gasProps1 = fluid.gasProps1();
    const gasProperties& gasProps2 = fluid.gasProps2();

    Field<scalar> drhodp1(mesh.nCells());
    Field<scalar> drhodp2(mesh.nCells());
    Field<scalar> drhodT1(mesh.nCells());
    Field<scalar> drhodT2(mesh.nCells());
    Field<scalar> H1(mesh.nCells());
    Field<scalar> H2(mesh.nCells());
    Field<scalar> cp1(mesh.nCells());
    Field<scalar> cp2(mesh.nCells());
    forAll(a1, celli)
    {
        drhodp1[celli] =  rho1[celli]*gasProps1.beta_T(p[celli], T1[celli]);
        drhodp2[celli] =  rho2[celli]*gasProps2.beta_T(p[celli], T2[celli]);
        drhodT1[celli] = -rho1[celli]*gasProps1.beta_p(p[celli], T1[celli]);
        drhodT2[celli] = -rho2[celli]*gasProps2.beta_p(p[celli], T2[celli]);
        H1[celli]  = e1[celli] + 0.5*magSqr(U1[celli]) + p[celli]/rho1[celli];
        H2[celli]  = e2[celli] + 0.5*magSqr(U2[celli]) + p[celli]/rho2[celli];
        cp1[celli] = gasProps1.Cp(p[celli], T1[celli]);
        cp2[celli] = gasProps2.Cp(p[celli], T2[celli]);
    }
    
    volScalarField::Internal dVdt = mesh.V()/dt;

    forAll(mesh.cells(), celli)
    {
        Eigen::FullPivLU<Eigen::Matrix<double,10,10>> D_inv{};

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

        D_inv.setThreshold(1e-300);
        D_inv.compute(J + pIntAdditionJ);// + (dragJ*mesh.V()[celli]));

        scalar rAlphaRho1  = /*(mesh.V()[celli]/dt.value())**/rezAlphaRho1[celli];
        scalar rAlphaRho2  = /*(mesh.V()[celli]/dt.value())**/rezAlphaRho2[celli];
        vector rAlphaRhoU1 = /*(mesh.V()[celli]/dt.value())**/rezAlphaRhoU1[celli];
        vector rAlphaRhoU2 = /*(mesh.V()[celli]/dt.value())**/rezAlphaRhoU2[celli];
        scalar rEpsilon1   = /*(mesh.V()[celli]/dt.value())**/rezEpsilon1[celli];
        scalar rEpsilon2   = /*(mesh.V()[celli]/dt.value())**/rezEpsilon2[celli];
    
        Eigen::Vector<double,10> b;
        b << rAlphaRho1, rAlphaRhoU1.x(), rAlphaRhoU1.y(), rAlphaRhoU1.z(), rEpsilon1,  rAlphaRho2, rAlphaRhoU2.x(), rAlphaRhoU2.y(), rAlphaRhoU2.z(), rEpsilon2;

        Eigen::Vector<double,10> starState = D_inv.solve(b);

        dp[celli]     = dt.value()*starState(0);
        dalpha[celli] = dt.value()*starState(1);
        dU1[celli]    = dt.value()*vector(starState(2), starState(3), starState(4));
        dU2[celli]    = dt.value()*vector(starState(5), starState(6), starState(7));
        dT1[celli]    = dt.value()*starState(8);
        dT2[celli]    = dt.value()*starState(9);
    }
}
