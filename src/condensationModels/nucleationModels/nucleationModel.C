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

#include "nucleationModel.H"
#include "twoFluid.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
    defineTypeNameAndDebug(nucleationModel, 0);
    defineRunTimeSelectionTable(nucleationModel, dictionary);
}
}

const Foam::dimensionSet Foam::TwoFluidFoam::nucleationModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::nucleationModel::nucleationModel
(
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const volScalarField& sigma,
    const volScalarField& Ts,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, twoFluid.name()),
            typeName,
            gasThermo.mesh().time().timeName(),
            gasThermo.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    gasThermo_(gasThermo),
    liquidThermo_(liquidThermo),
    sigma_(sigma),
    Ts_(Ts)
{}


Foam::TwoFluidFoam::nucleationModel::nucleationModel
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const volScalarField& sigma,
    const volScalarField& Ts,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, fluid.name()),
            typeName,
            gasThermo.mesh().time().timeName(),
            gasThermo.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    gasThermo_(gasThermo),
    liquidThermo_(liquidThermo),
    sigma_(sigma),
    Ts_(Ts)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::TwoFluidFoam::nucleationModel>
Foam::TwoFluidFoam::nucleationModel::New
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const volScalarField& sigma,
    const volScalarField& Ts
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting nucleationModel for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "nucleationModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, gasThermo, liquidThermo, sigma, Ts, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::nucleationModel::~nucleationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::nucleationModel::J(const volScalarField& rc) const
{
    const volScalarField sigma;

    const volScalarField& rho_g = gasThermo_.rho();
    const volScalarField& rho_l = liquidThermo_.rho();
    const volScalarField& T_g = gasThermo_.T();

    volScalarField J = sqrt(2*sigma/(pi*pow3(m1_)))*sqr(rho_g)/rhos_l*exp(-beta_*4*pi*sqr(rc)*sigma/(3*kB*T_g));

    return pos(Ts_ - T_g)*J;
}


bool Foam::TwoFluidFoam::nucleationModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
