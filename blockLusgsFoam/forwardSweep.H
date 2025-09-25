{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();
        const vectorField& pSf = mesh.Sf().boundaryField()[patchi];
        const scalarField& pMeshPhi = meshPhi.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            label ci = pFaceCells[facei];
            scalar ac = 0.5*lusgsOmega*(mag(U[ci] & pSf[facei]) + a[ci]*mag(pSf[facei]));
            D[ci] += ac;
        }
    }

    forAll(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rRho  = mesh.V()[cellI]/dt.value()*dRho[cellI];
        vector rRhoU = mesh.V()[cellI]/dt.value()*dRhoU[cellI];
        scalar rRhoE = mesh.V()[cellI]/dt.value()*dRhoE[cellI];

        forAll(cellFaces,i)
        {
            label faceI = cellFaces[i];                      // face index

            if (mesh.isInternalFace(faceI))
            {
                // Calculate D
                label own = owner[faceI];
                label nei = neighbour[faceI];

                scalar ac = 0.5*lusgsOmega
                    *(mag(U[cellI] & Sf[faceI]) + a[cellI]*magSf[faceI]);

                D[cellI] += ac;

                // we use the mesh ordering property: own < nei
                if (own != cellI)  // then data[own] belongs to lower triangle
                {
                     // Calculate L based on A+
                    scalar ac = 0.5*lusgsOmega*( 
                        mag(U[own] & Sf[faceI]) + a[own]*magSf[faceI]);
                    
                    scalar rho1  = rho[own]  + dRhoStar[own];
                    vector rhoU1 = rhoU[own] + dRhoUStar[own];
                    scalar rhoE1 = rhoE[own] + dRhoEStar[own];
                    
                    scalar e1 = rhoE1/rho1 - 0.5*magSqr(rhoU1)/sqr(rho1);
                    scalar p1 = gasProps.pEs(e1, rho1, p[own]);

                    scalar phi0 = U[own] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];

                    rRho  += ac*dRhoStar[own]  + 0.5*(rho1*phi1 - rho[own]*phi0); 
                    rRhoU += ac*dRhoUStar[own] + 0.5*((rhoU1*phi1 - rhoU[own]*phi0) + (p1 - p[own])*Sf[faceI]);
                    rRhoE += ac*dRhoEStar[own] + 0.5*(rhoE1*phi1 + p1*phi1 - (rhoE[own]*phi0 + p[own]*phi0)); 
                }
            }
        }

        dRhoStar[cellI]  = rRho  / D[cellI];
        dRhoUStar[cellI] = rRhoU / D[cellI];
        dRhoEStar[cellI] = rRhoE / D[cellI];
    }
}
