/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "totalLoss.H"
#include "fluidThermo.H"
#include "fvcVolumeIntegrate.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalLoss, 0);
    addToRunTimeSelectionTable(functionObject, totalLoss, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::totalLoss::totalLoss
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::totalLoss::read(const dictionary& dict)
{
    return fvMeshFunctionObject::read(dict);
}


bool Foam::functionObjects::totalLoss::execute()
{
    return true;
}


bool Foam::functionObjects::totalLoss::write()
{
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volScalarField& T1 = mesh_.lookupObject<volScalarField>("T.1");
    const volScalarField& T2 = mesh_.lookupObject<volScalarField>("T.2");

    const fluidThermo& thermo1 =
        lookupObject<fluidThermo>(IOobject::groupName("thermophysicalProperties", "1"));
    
    const fluidThermo& thermo2 =
        lookupObject<fluidThermo>(IOobject::groupName("thermophysicalProperties", "2"));

    const volScalarField h1 = thermo1.he() + p/thermo1.rho();
    const volScalarField h2 = thermo2.he() + p/thermo2.rho();

    const volScalarField& mDot = mesh_.lookupObject<volScalarField>("mDotCondensation");

    const volScalarField zRelax
    (
        IOobject
        (
            scopedName("zRelax"),
            time_.timeName(),
            mesh_
        ),
        mDot*(h1 - h2)*(1.0/T1 - 1.0/T2)*(T1 + T2)/2.0
    );

    Log << "    Relaxation losses field " << zRelax.name()
        << " to " << time_.timeName() << ", total losses: " << fvc::domainIntegrate(zRelax) << endl;

    zRelax.write();

    return true;
}


// ************************************************************************* //

