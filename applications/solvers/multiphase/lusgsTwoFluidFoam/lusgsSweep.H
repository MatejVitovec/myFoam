{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    const gasProperties& gasProps1 = fluid.gas1();
    const gasProperties& gasProps2 = fluid.gas2();

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

    Field<scalar> Dinv_p_alphaRho1 (mesh.nCells());
    Field<vector> Dinv_p_alphaRhoU1(mesh.nCells());
    Field<scalar> Dinv_p_alphaRhoE1(mesh.nCells());
    Field<scalar> Dinv_p_alphaRho2 (mesh.nCells());
    Field<vector> Dinv_p_alphaRhoU2(mesh.nCells());
    Field<scalar> Dinv_p_alphaRhoE2(mesh.nCells());

    Field<scalar> Dinv_alpha_alphaRho1 (mesh.nCells());
    Field<vector> Dinv_alpha_alphaRhoU1(mesh.nCells());
    Field<scalar> Dinv_alpha_alphaRhoE1(mesh.nCells());
    Field<scalar> Dinv_alpha_alphaRho2 (mesh.nCells());
    Field<vector> Dinv_alpha_alphaRhoU2(mesh.nCells());
    Field<scalar> Dinv_alpha_alphaRhoE2(mesh.nCells());

    Field<vector> Dinv_U1_alphaRho1 (mesh.nCells());
    Field<tensor> Dinv_U1_alphaRhoU1(mesh.nCells());
    Field<vector> Dinv_U1_alphaRhoE1(mesh.nCells());
    Field<vector> Dinv_U1_alphaRho2 (mesh.nCells());
    Field<tensor> Dinv_U1_alphaRhoU2(mesh.nCells());
    Field<vector> Dinv_U1_alphaRhoE2(mesh.nCells());

    Field<vector> Dinv_U2_alphaRho1 (mesh.nCells());
    Field<tensor> Dinv_U2_alphaRhoU1(mesh.nCells());
    Field<vector> Dinv_U2_alphaRhoE1(mesh.nCells());
    Field<vector> Dinv_U2_alphaRho2 (mesh.nCells());
    Field<tensor> Dinv_U2_alphaRhoU2(mesh.nCells());
    Field<vector> Dinv_U2_alphaRhoE2(mesh.nCells());
    
    Field<scalar> Dinv_T1_alphaRho1 (mesh.nCells());
    Field<vector> Dinv_T1_alphaRhoU1(mesh.nCells());
    Field<scalar> Dinv_T1_alphaRhoE1(mesh.nCells());
    Field<scalar> Dinv_T1_alphaRho2 (mesh.nCells());
    Field<vector> Dinv_T1_alphaRhoU2(mesh.nCells());
    Field<scalar> Dinv_T1_alphaRhoE2(mesh.nCells());

    Field<scalar> Dinv_T2_alphaRho1 (mesh.nCells());
    Field<vector> Dinv_T2_alphaRhoU1(mesh.nCells());
    Field<scalar> Dinv_T2_alphaRhoE1(mesh.nCells());
    Field<scalar> Dinv_T2_alphaRho2 (mesh.nCells());
    Field<vector> Dinv_T2_alphaRhoU2(mesh.nCells());
    Field<scalar> Dinv_T2_alphaRhoE2(mesh.nCells());
    
    #include "calculateDinv.H"

    Field<scalar> dpStar(mesh.nCells());
    Field<scalar> dalphaStar(mesh.nCells());
    Field<vector> dU1Star(mesh.nCells());
    Field<vector> dU2Star(mesh.nCells());
    Field<scalar> dT1Star(mesh.nCells());
    Field<scalar> dT2Star(mesh.nCells());
    
    #include "forwardSweep.H"
    #include "backwardSweep.H"
}
