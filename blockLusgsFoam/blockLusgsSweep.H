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

    //Field<scalar> D(mesh.V()/dt);
    /*if (LTS) {
        D *= dt * trDeltaTau();
    }*/

    Field<scalar> a(mesh.cells().size());
    forAll (mesh.cells(), celli)
    {
        a[celli] = gasProps.c(p[celli], T[celli]);
    }

    //#include "calculateDinv.H"

    Field<scalar> dRhoStar (dRho);
    Field<vector> dRhoUStar(dRhoU);
    Field<scalar> dRhoEStar(dRhoE);

    Field<scalar> dpStar(dp);
    Field<vector> dUStar(dU);
    Field<scalar> dTStar(dT);
    
    #include "blockForwardSweep.H"
    #include "blockBackwardSweep.H"
}
