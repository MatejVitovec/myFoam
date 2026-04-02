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

    forAll(owner, faceI)
    {

        /*scalar cm_own = 1.0 / Foam::sqrt((1.0 - alpha[owner[faceI]])/(rho1[owner[faceI]]*a1[owner[faceI]]*a1[owner[faceI]]) + alpha[owner[faceI]]/(rho2[owner[faceI]]*a2[owner[faceI]]*a2[owner[faceI]]));
        scalar cm_nei = 1.0 / Foam::sqrt((1.0 - alpha[neighbour[faceI]])/(rho1[neighbour[faceI]]*a1[neighbour[faceI]]*a1[neighbour[faceI]]) + alpha[neighbour[faceI]]/(rho2[neighbour[faceI]]*a2[neighbour[faceI]]*a2[neighbour[faceI]]));

        ac[faceI] = lusgsOmega*max(max(mag(U1[owner[faceI]]     & Sf[faceI]) + cm_own*magSf[faceI],
                                       mag(U1[neighbour[faceI]] & Sf[faceI]) + cm_nei*magSf[faceI]),
                                   max(mag(U2[owner[faceI]]     & Sf[faceI]) + cm_own*magSf[faceI],
                                       mag(U2[neighbour[faceI]] & Sf[faceI]) + cm_nei*magSf[faceI]));*/
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        forAll(mesh.boundary()[patchi], faceI)
        {
            /*ac[pFaceCells[faceI]] = lusgsOmega*max(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + a1[pFaceCells[faceI]]*magSf[faceI],
                                                   mag(U2[pFaceCells[faceI]] & Sf[faceI]) + a2[pFaceCells[faceI]]*magSf[faceI]);*/

            //ac[pFaceCells[faceI]] = lusgsOmega*(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + a1[pFaceCells[faceI]]*magSf[faceI]);

            /*scalar cm = 1.0 / Foam::sqrt((1.0 - alpha[pFaceCells[faceI]])/(rho1[pFaceCells[faceI]]*a1[pFaceCells[faceI]]*a1[pFaceCells[faceI]]) + alpha[pFaceCells[faceI]]/(rho2[pFaceCells[faceI]]*a2[pFaceCells[faceI]]*a2[pFaceCells[faceI]]));

            ac[pFaceCells[faceI]] = lusgsOmega*max(mag(U1[pFaceCells[faceI]] & Sf[faceI]) + cm*magSf[faceI],
                                                   mag(U2[pFaceCells[faceI]] & Sf[faceI]) + cm*magSf[faceI]);*/
        }
    }

    //std::vector<Eigen::FullPivLU<Eigen::Matrix<double,10,10>>> D_inv(mesh.nCells());
    std::vector<Eigen::Matrix<double,10,10>> J_(mesh.nCells());
    std::vector<Eigen::ColPivHouseholderQR<Eigen::Matrix<double,10,10>>> D_inv(mesh.nCells());

    //std::vector<scalar> alpha1zero(mesh.nCells());
    //std::vector<scalar> alpha2zero(mesh.nCells());

    //SVD
    /*std::vector<std::unique_ptr<Eigen::JacobiSVD<Eigen::Matrix<double,10,10>>>> D_inv(mesh.nCells());
    for(int celli = 0; celli < D_inv.size(); celli++)
    {
        D_inv[celli] = (std::make_unique<Eigen::JacobiSVD<Eigen::Matrix<double,10,10>>>());
    }*/

    
    #include "calculateDinv.H"

    volScalarField::Internal dpStar(IOobject("dpStar", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("dpStar", dimPressure, 0.0));
    volScalarField::Internal dalphaStar(IOobject("dalphaStar", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("dalphaStar", dimless, 0.0));
    volVectorField::Internal dU1Star(IOobject("dUStar.1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedVector("dU1", dimVelocity, vector(0, 0, 0)));
    volVectorField::Internal dU2Star(IOobject("dUStar.2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedVector("dU1", dimVelocity, vector(0, 0, 0)));
    volScalarField::Internal dT1Star(IOobject("dTStar.1", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("dT1Star", dimTemperature, 0.0));
    volScalarField::Internal dT2Star(IOobject("dTStar.2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("dT2Star", dimTemperature, 0.0));

    volScalarField pStar(IOobject("pStar", p.time().timeName(), p.db(), IOobject::NO_READ, IOobject::NO_WRITE), p);
    volScalarField alphaStar(IOobject("alphaStar", alpha.time().timeName(), alpha.db(), IOobject::NO_READ, IOobject::NO_WRITE), alpha);
    volVectorField U1Star(IOobject("UStar.1", U1.time().timeName(), U1.db(), IOobject::NO_READ, IOobject::NO_WRITE), U1);
    volVectorField U2Star(IOobject("UStar.2", U2.time().timeName(), U2.db(), IOobject::NO_READ, IOobject::NO_WRITE), U2);
    volScalarField T1Star(IOobject("TStar.1", T1.time().timeName(), T1.db(), IOobject::NO_READ, IOobject::NO_WRITE), T1);
    volScalarField T2Star(IOobject("TStar.2", T2.time().timeName(), T2.db(), IOobject::NO_READ, IOobject::NO_WRITE), T2);

    for (int intIter = 0; intIter < lusgsIntIters; intIter++)
    {
        #include "forwardSweep.H"
        //#include "boundarySweep.H"
        #include "backwardSweep.H"
    }
}
