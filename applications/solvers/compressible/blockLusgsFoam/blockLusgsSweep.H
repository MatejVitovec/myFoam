{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    //Field<scalar> drhodp(mesh.nCells());
    //Field<scalar> drhodT(mesh.nCells());
    //Field<scalar> H(mesh.nCells());
    //Field<scalar> cp(mesh.nCells());

    //Field<scalar> subA(mesh.nCells());
    //Field<scalar> subB(mesh.nCells());
    //Field<scalar> subC(mesh.nCells());

    volScalarField::Internal H = e + 0.5*magSqr(U) + p/rho;
    volScalarField::Internal drhodp = rho/p; //TODO netreba kopirovat pouze pro jednotky
    volScalarField::Internal drhodT = rho/T;
    volScalarField::Internal cp = thermo.Cp();

    forAll(H, celli)
    {
        drhodp[celli] =  rho[celli]*gasProps.beta_T(p[celli], T[celli]);
        drhodT[celli] = -rho[celli]*gasProps.beta_p(p[celli], T[celli]);
        //H[celli]  = e[celli] + 0.5*magSqr(U[celli]) + p[celli]/rho[celli];
        //cp[celli] = gasProps.Cp(p[celli], T[celli]);

        //subA =  drhodp[celli]*rho[celli]*rho[celli]*cp[celli] - drhodT[celli]*drhodT[celli]*T[celli];
        //subB = -drhodp[celli]*rho[celli]*(H[celli] - magSqr(U[celli])) - drhodT[celli]*T[celli];
        //subC =  drhodT[celli]*rho[celli]*(H[celli] - magSqr(U[celli])) + cp[celli]*rho[celli]*rho[celli];
    }

    volScalarField::Internal subA =  drhodp*rho*rho*cp - drhodT*drhodT*T;
    volScalarField::Internal subB = -drhodp*rho*(H - magSqr(U)().internalField()) - drhodT*T;
    volScalarField::Internal subC =  drhodT*rho*(H - magSqr(U)().internalField()) + cp*rho*rho;

    volScalarField::Internal D = mesh.V()/dt;

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

    Field<scalar> dpStar(mesh.nCells());
    Field<vector> dUStar(mesh.nCells());
    Field<scalar> dTStar(mesh.nCells());

    for (int intIter = 0; intIter < lusgsIntIters; intIter++)
    {
        #include "blockForwardSweep.H"
        #include "blockBackwardSweep.H"
    }
}
