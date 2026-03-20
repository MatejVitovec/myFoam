/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "virtualMassModel.H"
#include "phasePair.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
    defineTypeNameAndDebug(virtualMassModel, 0);
    defineRunTimeSelectionTable(virtualMassModel, dictionary);
}
}

const Foam::dimensionSet Foam::virtualMassModel::dimK(dimDensity);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::virtualMassModel::virtualMassModel
(
    const twoFluid& fluid,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, pair.name()),
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


Foam::virtualMassModel::virtualMassModel
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
            //IOobject::groupName(typeName, pair.name()),
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

Foam::autoPtr<Foam::virtualMassModel>
Foam::virtualMassModel::New
(
    const twoFluid& fluid
)
{
    const dictionary& dict = fluid.subDict("virtualMass");

    Info<< "Selecting virtualMassModel for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "virtualMassModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, fluid, true);
}


Foam::autoPtr<Foam::virtualMassModel>
Foam::virtualMassModel::New
(
    const dictionary& dict,
    const twoFluid& fluid
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting virtualMassModel for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "virtualMassModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, fluid, true);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::virtualMassModel::Ki() const
{
    return Cvm()*fluid_.thermo1().rho();
}


Foam::tmp<Foam::volScalarField> Foam::virtualMassModel::K() const
{
    return fluid_.alpha()*Ki();
}


Foam::tmp<Foam::surfaceScalarField> Foam::virtualMassModel::Kf() const
{
    return fvc::interpolate(fluid_.alpha())*fvc::interpolate(Ki());
}


bool Foam::virtualMassModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
