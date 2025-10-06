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

    
    /*Field<scalar> pCond_rho_p(mesh.nCells());
    Field<vector> pCond_rho_U(mesh.nCells());
    Field<scalar> pCond_rho_T(mesh.nCells());

    Field<vector> pCond_rhoU_p(mesh.nCells());
    Field<tensor> pCond_rhoU_U(mesh.nCells());
    Field<vector> pCond_rhoU_T(mesh.nCells());
    
    Field<scalar> pCond_rhoE_p(mesh.nCells());
    Field<vector> pCond_rhoE_U(mesh.nCells());
    Field<scalar> pCond_rhoE_T(mesh.nCells());*/

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
