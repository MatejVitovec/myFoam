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

#include "dragModel.H"
#include "twoFluid.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, dictionary);
}
}

const Foam::dimensionSet Foam::TwoFluidFoam::dragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModel::dragModel
(
    const twoFluid& fluid,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, twoFluid.name()),
            typeName,
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    fluid_(fluid)
{}


Foam::TwoFluidFoam::dragModel::dragModel
(
    const dictionary& dict,
    const twoFluid& fluid,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, fluid.name()),
            typeName,
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    fluid_(fluid)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
Foam::autoPtr<Foam::TwoFluidFoam::dragModel>
Foam::TwoFluidFoam::dragModel::New
(
    const twoFluid& fluid
)
{
    const dictionary& dict = fluid.subDict("drag");

    const word modelType(dict.get<word>("type"));

    Info<< "Selecting dragModel for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "dragModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, fluid, true);
}


Foam::autoPtr<Foam::TwoFluidFoam::dragModel>
Foam::TwoFluidFoam::dragModel::New
(
    const dictionary& dict,
    const twoFluid& fluid
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting dragModel for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "dragModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, fluid, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModel::~dragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::dragModel::Ki(const volScalarField& d) const
{
    return (0.75*CdRe()*fluid_.thermo1().rho()*mag(fluid_.U1() - fluid_.U2()))/d;
}


Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::dragModel::K(const volScalarField& d) const
{
    return Ki(d)*fluid_.alpha();
}


bool Foam::TwoFluidFoam::dragModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
