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

#include "fvMesh.H"
#include "volFields.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> applyFunction1
(
    const std::function<scalar(scalar)> f,
    const volScalarField& src,
    const word name,
    const dimensionSet dim
)
{
    const fvMesh& mesh = src.mesh();
    tmp<volScalarField> tres
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dim
            )
        );
    volScalarField& r = tres.ref();

    forAll(r, celli)
    {
        r[celli] = f(src[celli]);
    }
    
    forAll(src.boundaryField(), patchi)
    {
        scalarField& rp = r.boundaryFieldRef()[patchi];
        const scalarField& sp = src.boundaryField()[patchi];
        forAll(rp, facei)
        {
            rp[facei] = f(sp[facei]);
        }
    }
    
    return tres;
}


tmp<volScalarField> applyFunction2
(
    const std::function<scalar(scalar,scalar)> f,
    const volScalarField& src1,
    const volScalarField& src2,
    const word name,
    const dimensionSet dim
)
{
    const fvMesh& mesh = src1.mesh();
    tmp<volScalarField> tres
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dim
            )
        );
    volScalarField& r = tres.ref();

    forAll(r, celli)
    {
        r[celli] = f(src1[celli], src2[celli]);
    }
    
    forAll(src1.boundaryField(), patchi)
    {
        scalarField& rp = r.boundaryFieldRef()[patchi];
        const scalarField& s1p = src1.boundaryField()[patchi];
        const scalarField& s2p = src2.boundaryField()[patchi];
        forAll(rp, facei)
        {
            rp[facei] = f(s1p[facei], s2p[facei]);
        }
    }
    
    return tres;
}

} // End namespace Foam

// ************************************************************************* //