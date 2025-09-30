{
    volScalarField::Internal D = mesh.V()/dt;

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();


    forAll(owner, faceI)
    {
        D[owner[faceI]]     += 0.5*lusgsOmega*(mag(U[owner[faceI]]     & Sf[faceI]) + a[owner[faceI]]    *magSf[faceI]);
        D[neighbour[faceI]] += 0.5*lusgsOmega*(mag(U[neighbour[faceI]] & Sf[faceI]) + a[neighbour[faceI]]*magSf[faceI]);
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], faceI)
        {
            D[pFaceCells[faceI]] += 0.5*lusgsOmega*(mag(U[pFaceCells[faceI]] & Sf[faceI]) + a[pFaceCells[faceI]]*magSf[faceI]);
        }
    }


    //Field<scalar> D(mesh.V()/dt);
    /*if (LTS) {
        D *= dt * trDeltaTau();
    }*/

    Field<scalar> a(mesh.cells().size());
    forAll (mesh.cells(), celli)
    {
        a[celli] = gasProps.c(p[celli], T[celli]);
    }

    Field<scalar> dRhoStar (dRho);
    Field<vector> dRhoUStar(dRhoU);
    Field<scalar> dRhoEStar(dRhoE);
    
    #include "blockForwardSweep.H"
    #include "blockBackwardSweep.H"
}
