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

#include "SchillerNaumann.H"
#include "twoFluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
namespace dragModels
{
    defineTypeNameAndDebug(SchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, SchillerNaumann, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModels::SchillerNaumann::SchillerNaumann
(
    const dictionary& dict,
    const twoFluid& fluid,
    const bool registerObject
)
:
    dragModel(dict, fluid, registerObject),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModels::SchillerNaumann::~SchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::dragModels::SchillerNaumann::CdRe(const volScalarField& d) const
{
    tmp<volScalarField> Re = mag(fluid_.U1() - fluid_.U2())*d/fluid_.thermo1().nu();

    return neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687)) + pos0(Re - 1000)*0.44*max(Re, residualRe_);
}


// ************************************************************************* //
