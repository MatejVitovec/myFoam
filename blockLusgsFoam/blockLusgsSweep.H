{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    volScalarField::Internal D = mesh.V()/dt;

    Field<scalar> Dinv_p_rho (mesh.nCells());
    Field<vector> Dinv_p_rhoU(mesh.nCells());
    Field<scalar> Dinv_p_rhoE(mesh.nCells());

    Field<vector> Dinv_U_rho (mesh.nCells());
    Field<tensor> Dinv_U_rhoU(mesh.nCells());
    Field<vector> Dinv_U_rhoE(mesh.nCells());
    
    Field<scalar> Dinv_T_rho (mesh.nCells());
    Field<vector> Dinv_T_rhoU(mesh.nCells());
    Field<scalar> Dinv_T_rhoE(mesh.nCells());

    
    /*Field<scalar> a(mesh.cells().size());
    forAll (mesh.cells(), celli)
    {
        a[celli] = gasProps.c(p[celli], T[celli]);
    }*/

    #include "calculateDinv.H"


    Field<scalar> dpStar(mesh.nCells());
    Field<vector> dUStar(mesh.nCells());
    Field<scalar> dTStar(mesh.nCells());

    #include "blockForwardSweep.H"
    #include "blockBackwardSweep.H"
}
