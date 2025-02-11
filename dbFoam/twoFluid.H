/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvMesh.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace twoFluid
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void primitiveFromConservative
(
    volScalarField& p,
    volScalarField& alpha,
    volVectorField& U1,
    volVectorField& U2,
    volScalarField& T1,
    volScalarField& T2,
    const volScalarField& alphaRho1,
    const volScalarField& alphaRho2,
    const volVectorField& alphaRhoU1,
    const volVectorField& alphaRhoU2,
    const volScalarField& alphaRhoE1,
    const volScalarField& alphaRhoE2,
    const psiThermo& thermo1,
    const psiThermo& thermo2,
    const fvMesh& mesh
)
{
    const auto& cells = mesh.C();

    forAll(cells, celli)
    {
        //TODO
    }
}

void blend
(
    volScalarField& p,
    volScalarField& alpha,
    volVectorField& U1,
    volVectorField& U2,
    volScalarField& T1,
    volScalarField& T2,
)
{
    //TODO
}


void updateIntPressure
(
    volScalarField& pInt,
    const volScalarField& alpha,
    const volVectorField& U1,
    const volVectorField& U2,
    const psiThermo& thermo1,
    const psiThermo& thermo2
)
{
    const scalar sigma = 2.0;
    const scalar epsilonP = 0.01;

    const auto& rho1 = thermo1.rho();
    const auto& rho2 = thermo2.rho();

    pInt = sigma*((alpha*rho1*(1.0 - alpha)*rho2)/
        (alpha*rho2 + (1.0 - alpha)*rho1))*magSqr(U2 - U1);
}




namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
div
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
)
{
    return tmp<GeometricField<Type, fvPatchField, volMesh>>
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            "div("+ssf.name()+')',
            surfaceIntegrate(ssf_pos, ssf_neg)
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf.dimensions()/dimVol, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    surfaceIntegrate(vf.primitiveFieldRef(), ssf_pos, ssf_neg);
    vf.correctBoundaryConditions();

    return tvf;
}


template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
)
{
    const fvMesh& mesh = ssf.mesh();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const Field<Type>& issf_pos = ssf_pos;
    const Field<Type>& issf_neg = ssf_neg;

    forAll(owner, facei)
    {
        ivf[owner[facei]] += issf_pos[facei];
        ivf[neighbour[facei]] -= issf_neg[facei];
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundary()[patchi].faceCells();

        const fvsPatchField<Type>& pssf_pos = ssf_pos.boundaryField()[patchi];

        forAll(mesh.boundary()[patchi], facei)
        {
            ivf[pFaceCells[facei]] += pssf_pos[facei];
        }
    }

    ivf /= mesh.Vsc();
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoFluid

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
