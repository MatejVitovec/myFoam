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
//TODO correct include

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TwoFluidFoam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void surfaceIntegrate
(
    Field<Type>& ivf,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
)
{
    const fvMesh& mesh = ssf_pos.mesh();

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


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
surfaceIntegrate
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_pos,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf_neg
)
{
    const fvMesh& mesh = ssf_pos.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh>> tvf
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "surfaceIntegrate("+ssf_pos.name()+')',
                ssf_pos.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>(ssf_pos.dimensions()/dimVol, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );
    GeometricField<Type, fvPatchField, volMesh>& vf = tvf.ref();

    surfaceIntegrate(vf.primitiveFieldRef(), ssf_pos, ssf_neg);
    vf.correctBoundaryConditions();

    return tvf;
}


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
            "div("+ssf_pos.name()+')',
            surfaceIntegrate(ssf_pos, ssf_neg)
        )
    );
}

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
div
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf_pos,
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf_neg
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh>> Div(fvc::div(tssf_pos(), tssf_neg()));
    tssf_pos.clear();
    tssf_neg.clear();
    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TwoFluidFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
