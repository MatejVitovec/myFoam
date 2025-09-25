{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    
    forAllReverse(mesh.cells(), cellI)
    {
        const labelList& cellFaces = mesh.cells()[cellI];    // list of faces enclosing cellI

        scalar rRho  = 0;
        vector rRhoU = vector(0, 0, 0);
        scalar rRhoE = 0;
            
        forAll(cellFaces, i)
        {
            label faceI = cellFaces[i];                      // face index
            if (mesh.isInternalFace(faceI))
            {
                label nei = neighbour[faceI];
                label own = owner[faceI];

                // we use the mesh ordering property: own < nei
                if (nei != cellI)  // then data[nei] belongs to upper triangle
                {
                    scalar ac = 0.5*lusgsOmega*( 
                        mag(U[nei] & Sf[faceI]) + a[nei]*magSf[faceI]);

                    scalar rho1  = rho[nei]  + dRho[nei];
                    vector rhoU1 = rhoU[nei] + dRhoU[nei];
                    scalar rhoE1 = rhoE[nei] + dRhoE[nei];
                    scalar e1 = rhoE1/rho1 - 0.5*magSqr(rhoU1)/sqr(rho1);
                    scalar p1 = gasProps.pEs(e1, rho1, p[nei]);

                    scalar phi0 = U[nei] & Sf[faceI];
                    scalar phi1 = (rhoU1/rho1) & Sf[faceI];
                    
                    rRho  += ac*dRho[nei]  - 0.5*(rho1*phi1 - rho[nei]*phi0); 
                    rRhoU += ac*dRhoU[nei] - 0.5*((rhoU1*phi1 - rhoU[nei]*phi0) + (p1 - p[nei])*Sf[faceI]);
                    rRhoE += ac*dRhoE[nei] - 0.5 * ((rhoE1*phi1 + p1*phi1) - (rhoE[nei]*phi0 + p[nei]*phi0)); 
                }
            }
        }

        dRho[cellI]  = dRhoStar[cellI]  + rRho  / D[cellI];
        dRhoU[cellI] = dRhoUStar[cellI] + rRhoU / D[cellI];
        dRhoE[cellI] = dRhoEStar[cellI] + rRhoE / D[cellI];
    }
}
