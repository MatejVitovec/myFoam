{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    volScalarField::Internal Dinv_p_p;
    volVectorField::Internal Dinv_p_U;
    volScalarField::Internal Dinv_p_T;

    volVectorField::Internal Dinv_U_p;
    volTensorField::Internal Dinv_U_U;
    volVectorField::Internal Dinv_U_T;

    volScalarField::Internal Dinv_T_p;
    volVectorField::Internal Dinv_T_U;
    volScalarField::Internal Dinv_T_T;

    forAll(mesh.cells(), cellI)
    {
        

        D[owner[faceI]]     += 0.5*lusgsOmega*(mag(U[owner[faceI]]     & Sf[faceI]) + a[owner[faceI]]    *magSf[faceI]);
        D[neighbour[faceI]] += 0.5*lusgsOmega*(mag(U[neighbour[faceI]] & Sf[faceI]) + a[neighbour[faceI]]*magSf[faceI]);
    }
}
