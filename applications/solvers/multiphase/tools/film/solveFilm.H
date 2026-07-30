{
    label wallPatchID = mesh.boundaryMesh().findPatchID("wall");

    if (wallPatchID < 0)
    {
        FatalErrorInFunction
            << "Patch wall not found"
            << abort(FatalError);
    }

    tmp<volSymmTensorField> tTau = turbulence1->devRhoReff();

    const scalarField& pAlphaRhoPhi2 = twoFluidFlux.alphaRhoFlux2_pos().boundaryField()[wallPatchID];

    const vectorField Up = U1.boundaryField()[wallPatchID].patchInternalField();
    const fvPatchScalarField& rhop = rho1.boundaryField()[wallPatchID];

    const fvPatch& patch = mesh.boundary()[wallPatchID];
    const vectorField pnf = patch.nf();

    dimensionedScalar rhoLiquid
    (
        "rhoLiquid",
        dimDensity,
        1000.0
    );

    dimensionedScalar muLiquid
    (
        "muLiquid",
        dimDynamicViscosity,
        1e-3
    );

    areaVectorField Upa
    (
        IOobject
        (
            "Ua",
            runTime.timeName(),
            aMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedVector(U1.dimensions(), Zero)
    );

    Upa.primitiveFieldRef() = Up;
    Upa.correctBoundaryConditions();

    areaScalarField Sfilm
    (
        IOobject
        (
            "Sfilm",
            runTime.timeName(),
            aMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    );

    areaScalarField tauG
    (
        IOobject
        (
            "tauG",
            runTime.timeName(),
            aMesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    );

    const scalar Cf = 0.01;

    forAll(h, faceI)
    {
        scalar Umag = mag(Up[faceI]);

        tauG[faceI] = 0.5*rhop[faceI]*sqr(Umag)*Cf;

        Sfilm[faceI] = pAlphaRhoPhi2[faceI]/rhoLiquid.value();
    }

    areaVectorField n = aMesh.faceAreaNormals();

    n /= mag(n);

    areaVectorField Ut = Upa - ((Upa & n)*n);

    areaVectorField tauWallCalc = (Ut/(mag(Ut) + dimensionedScalar("smallSpeed", dimVelocity, VSMALL)))*tauG;

    //Uf = (tauWallCalc*h)/(2.0*muLiquid);
    
    areaVectorField tauWall
    (
        IOobject
        (
            "tauWall",
            runTime.timeName(),
            aMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedVector(dimPressure, Zero)
    );

    vectorField traction = -(tTau().boundaryField()[wallPatchID] & pnf);
    tauWall.primitiveFieldRef() = (traction - pnf*(pnf & traction)); // pouze tečná složka = wall shear stress

    //areaVectorField UfAlt = (tauWall*h)/(2.0*muLiquid);

    Uf = (tauWall*h)/(2.0*muLiquid);

    /*forAll(Uf, faceI)
    {
        vector n = aMesh.faceAreaNormals()[faceI];

        n /= mag(n);

        vector Ut = Up[faceI] - (Up[faceI] & n)*n;

        const scalar magUt = mag(Ut);

        if (magUt > SMALL)
        {
            Uf[faceI] = ((Ut/magUt)*tauG[faceI]*h[faceI])/(2.0*muLiquid.value());
        }
        else
        {
            Uf[faceI] = Zero;
        }
    }*/

    edgeScalarField phiF
    (
        IOobject
        (
            "phiF",
            runTime.timeName(),
            aMesh
        ),
        fac::interpolate(Uf) & aMesh.Le()
    );

    faScalarMatrix hEqn
    (
        fam::ddt(h)
      + fam::div(phiF, h)
        ==
        Sfilm
    );

    hEqn.relax();

    hEqn.solve();

    h.max(0.0);


    volScalarField hWall
    (
        IOobject
        (
            "hWall",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimLength, Zero)
    );

    volVectorField uWall
    (
        IOobject
        (
            "uWall",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimLength/dimTime, vector::zero)
    );

    hWall.boundaryFieldRef()[wallPatchID] == h;
    uWall.boundaryFieldRef()[wallPatchID] == Uf;

    if (mesh.time().outputTime())
    {
        hWall.write();
        uWall.write();
    }

    Info << "Film solved" << endl;
}
