{
    Field<scalar> D(mesh.V()/dt);
    /*if (LTS) {
        D *= dt * trDeltaTau();
    }*/

    Field<scalar> a(mesh.cells().size());
    forAll (mesh.cells(), celli)
    {
        a[celli] = gasProps.c(p[celli], T[celli]);
    }

    Field<scalar> dRhoStar(dRho);
    Field<vector> dRhoUStar(dRhoU);
    Field<scalar> dRhoEStar(dRhoE);
    
    #include "blockForwardSweep.H"
    #include "blockBackwardSweep.H"
}
