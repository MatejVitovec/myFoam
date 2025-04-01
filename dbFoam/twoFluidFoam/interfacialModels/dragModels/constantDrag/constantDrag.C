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

#include "constantDrag.H"
#include "twoFluid.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantDrag, 0);
    defineRunTimeSelectionTable(constantDrag, dictionary);
}

const Foam::dimensionSet Foam::constantDrag::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantDrag::constantDrag
(
    const twoFluid& fluid,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            IOobject::groupName(typeName, twoFluid.name()),
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    fluid_(fluid)
{}


Foam::constantDrag::constantDrag
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
            IOobject::groupName(typeName, fluid.name()),
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

Foam::autoPtr<Foam::constantDrag>
Foam::constantDrag::New
(
    const dictionary& dict,
    const twoFluid& fluid
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting constantDrag for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "constantDrag",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, pair, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantDrag::~constantDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::constantDrag::Ki() const
{
    const dimensionedScalar d("dropletDiameter", dimLength, 1e-5);

    return
        0.75
       *CdRe()
       //*swarmCorrection_->Cs()
       *fluid_.thermo1().rho()
       //*fluid_.thermo1().nu()
       ///sqr(pair_.dispersed().d());
       /d;
}


Foam::tmp<Foam::volScalarField> Foam::constantDrag::K() const
{
    return (1.0 - fluid_.alpha())*Ki();
}


bool Foam::constantDrag::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
