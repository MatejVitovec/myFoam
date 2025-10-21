{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    Field<scalar> Dinv_p_rho (mesh.nCells());
    Field<vector> Dinv_p_rhoU(mesh.nCells());
    Field<scalar> Dinv_p_rhoE(mesh.nCells());

    Field<vector> Dinv_U_rho (mesh.nCells());
    Field<tensor> Dinv_U_rhoU(mesh.nCells());
    Field<vector> Dinv_U_rhoE(mesh.nCells());
    
    Field<scalar> Dinv_T_rho (mesh.nCells());
    Field<vector> Dinv_T_rhoU(mesh.nCells());
    Field<scalar> Dinv_T_rhoE(mesh.nCells());

    Field<scalar> drhodp(mesh.nCells());
    Field<scalar> drhodT(mesh.nCells());
    Field<scalar> H(mesh.nCells());
    Field<scalar> cp(mesh.nCells());
    forAll(a, celli)
    {
        drhodp[celli] =  rho[celli]*gasProps.beta_T(p[celli], T[celli]);
        drhodT[celli] = -rho[celli]*gasProps.beta_p(p[celli], T[celli]);
        H[celli]  = e[celli] + 0.5*magSqr(U[celli]) + p[celli]/rho[celli];
        cp[celli] = gasProps.Cp(p[celli], T[celli]);
    }

    #include "calculateDinv.H"

    Field<scalar> dpStar(mesh.nCells());
    Field<vector> dUStar(mesh.nCells());
    Field<scalar> dTStar(mesh.nCells());

    for (int intIter = 0; intIter < lusgsIntIters; intIter++)
    {
        #include "blockForwardSweep.H"
        #include "blockBackwardSweep.H"
    }
}
