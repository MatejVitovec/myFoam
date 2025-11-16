/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "growthModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{
    defineTypeNameAndDebug(growthModel, 0);
    defineRunTimeSelectionTable(growthModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WetSteam::growthModel::growthModel
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur
)
:
    gasThermo_(gasThermo),
    liquidThermo_(liquidThermo),
    saturation_(satur),
    rMin_
    (
        "rMin",
        dimLength,
        dict.lookupOrDefault<scalar>("rMin", 1e-9)
    ),
    rDot_
    (
        IOobject
        (
            "rDot",
            gasThermo_.p().mesh().time().timeName(),
            gasThermo_.p().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        gasThermo_.p().mesh(),
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::WetSteam::growthModel>
Foam::WetSteam::growthModel::New
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur
)
{
    word growthModelName(dict.lookup("growthModel"));

    Info<< "Selecting growthModel "
        << growthModelName << endl;

    auto* ctorPtr = dictionaryConstructorTable(growthModelName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "growthModel",
            growthModelName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, gasThermo, liquidThermo, satur);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::WetSteam::growthModel::~growthModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar Foam::WetSteam::growthModel::rMin() const
{
    return rMin_;
}


const Foam::volScalarField& Foam::WetSteam::growthModel::rDot() const
{
    return rDot_;
}

// ************************************************************************* //
