{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    const gasProperties& gasPropsv = fluid.gasProps1();
    const gasProperties& gasPropsl = fluid.gasProps2();

    const volScalarField& rhol = fluid.thermo2().rho();
    const volScalarField& rhov = fluid.thermo1().rho();

    const volScalarField& el = fluid.thermo2().he();
    const volScalarField& ev = fluid.thermo1().he();

    const volScalarField E = w*el + (1.0 - w)*ev + 0.5*magSqr(U);

    // prirazeni pouze pro jednotky TODO
    volScalarField::Internal drhodp = rho/p;;
    volScalarField::Internal drhodT = rho/T;
    volScalarField::Internal drhodw = rho;
    volScalarField::Internal dEdp   = E/p;
    volScalarField::Internal dEdT   = E/T;
    volScalarField::Internal dEdw   = E;
    forAll(rho, celli)
    {
        const scalar beta_Tl = gasPropsl.beta_T(p[celli], T[celli]);
        const scalar beta_Tv = gasPropsv.beta_T(p[celli], T[celli]);
        const scalar beta_pl = gasPropsl.beta_p(p[celli], T[celli]);
        const scalar beta_pv = gasPropsv.beta_p(p[celli], T[celli]);
        drhodp[celli] = ((w[celli]/sqr(rhol[celli]))*(rho2[celli]*beta_Tl)
               + ((1.0 - w[celli])/sqr(rhov[celli]))*(rho1[celli]*beta_Tv))
               /sqr(rho[celli]);
        drhodT[celli] = ((w[celli]/sqr(rhol[celli]))*(-rho2[celli]*beta_pl)
               + ((1.0 - w[celli])/sqr(rhov[celli]))*(-rho1[celli]*beta_pv))
               /sqr(rho[celli]);
        drhodw[celli] = (1.0/rhov[celli] - 1.0/rhol[celli])/sqr(rho[celli]);

        const scalar deldp = p[celli]/rhol[celli]*beta_pl - T[celli]/rhol[celli]*beta_Tl;
        const scalar devdp = p[celli]/rhov[celli]*beta_pv - T[celli]/rhov[celli]*beta_Tv;
        dEdp[celli] = w[celli]*deldp + (1.0 - w[celli]*devdp);

        const scalar deldT = gasPropsl.Cp(p[celli], T[celli]) - p[celli]/rhol[celli]*beta_pl;
        const scalar devdT = gasPropsv.Cp(p[celli], T[celli]) - p[celli]/rhov[celli]*beta_pv;
        dEdT[celli] = w[celli]*deldT + (1.0 - w[celli]*devdT);

        dEdw[celli] = el[celli] - ev[celli];
    }

    volScalarField::Internal invD = dEdT*drhodp*rho - dEdp*drhodT*rho;
    volScalarField::Internal invA = rho*dEdT + E*drhodT - drhodT*magSqr(U)().internalField() - dEdw*drhodT*w + dEdT*drhodw*w;
    volScalarField::Internal invB = rho*dEdp - E*drhodp - drhodp*magSqr(U)().internalField() - dEdw*drhodp*w + dEdp*drhodw*w;
    volScalarField::Internal invC =  dEdw*drhodT - dEdT*drhodw;
    volScalarField::Internal invE = -dEdw*drhodp + dEdp*drhodw;

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
    Field<scalar> dwStar(mesh.nCells());

    for (int intIter = 0; intIter < lusgsIntIters; intIter++)
    {
        #include "forwardSweep.H"
        #include "backwardSweep.H"
    }
}
