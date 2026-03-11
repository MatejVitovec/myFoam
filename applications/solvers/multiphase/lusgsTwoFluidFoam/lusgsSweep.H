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


    surfaceScalarField ac
    (
        IOobject
        (
            "ac",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("ac", dimLength/dimTime, 0.0)
    );

    forAll(owner, faceI)
    {
        ac[faceI] = lusgsOmega*max(max(mag(U1[owner[faceI]]     & Sf[faceI]) + a1[owner[faceI]]    *magSf[faceI],
                                       mag(U1[neighbour[faceI]] & Sf[faceI]) + a1[neighbour[faceI]]*magSf[faceI]),
                                   max(mag(U2[owner[faceI]]     & Sf[faceI]) + a2[owner[faceI]]    *magSf[faceI],
                                       mag(U2[neighbour[faceI]] & Sf[faceI]) + a2[neighbour[faceI]]*magSf[faceI]));
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], faceI)
        {
            ac[pFaceCells[faceI]] = lusgsOmega*max(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + a1[pFaceCells[faceI]]*magSf[faceI],
                                                   mag(U2[pFaceCells[faceI]] & Sf[faceI]) + a2[pFaceCells[faceI]]*magSf[faceI]);
        }
    }

    //std::vector<Eigen::FullPivLU<Eigen::Matrix<double,10,10>>> D_inv(mesh.nCells());
    std::vector<Eigen::ColPivHouseholderQR<Eigen::Matrix<double,10,10>>> D_inv(mesh.nCells()); 
    
    #include "calculateDinv.H"

    Field<scalar> dpStar(mesh.nCells());
    Field<scalar> dalphaStar(mesh.nCells());
    Field<vector> dU1Star(mesh.nCells());
    Field<vector> dU2Star(mesh.nCells());
    Field<scalar> dT1Star(mesh.nCells());
    Field<scalar> dT2Star(mesh.nCells());

    for (int intIter = 0; intIter < lusgsIntIters; intIter++)
    {
        #include "forwardSweep.H"
        #include "backwardSweep.H"
    }
}
