/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "constantDrag.H"
#include "twoFluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
namespace dragModels
{
    defineTypeNameAndDebug(constantDrag, 0);
    addToRunTimeSelectionTable(dragModel, constantDrag, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModels::constantDrag::constantDrag
(
    const dictionary& dict,
    const twoFluid& fluid,
    const bool registerObject
)
:
    dragModel(dict, fluid, registerObject),
    dragCoeff_("dragCoeff", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModels::constantDrag::~constantDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::dragModels::constantDrag::CdRe() const
{
    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "constantDragCoeff",
                    fluid_.mesh().time().timeName(),
                    fluid_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fluid_.mesh(),
                dragCoeff_
            )
        );
}


// ************************************************************************* //
