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

#include "viscousDrag.H"
#include "twoFluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
    defineTypeNameAndDebug(viscousDrag, 0);
    addToRunTimeSelectionTable(dragModel, viscousDrag, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::viscousDrag::viscousDrag
(
    const dictionary& dict,
    const twoFluid& fluid,
    const volScalarField& dropletDiameter,
    const bool registerObject
)
:
    dragModel(dict, fluid, dropletDiameter, registerObject)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::viscousDrag::~viscousDrag()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::viscousDrag::Cc(const volScalarField& d) const
{
    //const objectRegistry& db = fluid_.mesh();
    //const volScalarField& Kn = db.lookupObject<volScalarField>("Kn");

    const dimensionedScalar Rg("SpecificGasConstant", dimEnergy/dimMass/dimTemperature, 461.68449);
    volScalarField Kn = 1.5*fluid_.thermo1().mu()*sqrt(Rg*fluid_.T1())/(d*fluid_.p());

    //volScalarField Kn = 1.5*fluid_.thermo1().mu()*sqrt((Foam::constant::thermodynamic::RR/fluid_.thermo1().W())*fluid_.T1())/(d*fluid_.p());

    return 1.0 + 2*Kn*(1.257 + 0.4*exp(-1.1/(2.0*Kn + dimensionedScalar("dimlessNearZero", dimless, SMALL))));
}


Foam::scalar Foam::TwoFluidFoam::viscousDrag::dKidp(const label celli) const
{
    return 0.0;
}

Foam::scalar Foam::TwoFluidFoam::viscousDrag::dKidalpha(const label celli) const
{
    return 0.0;
}

Foam::vector Foam::TwoFluidFoam::viscousDrag::dKidU1(const label celli) const
{
    return vector::zero;
}

Foam::vector Foam::TwoFluidFoam::viscousDrag::dKidU2(const label celli) const
{
    return vector::zero;
}

Foam::scalar Foam::TwoFluidFoam::viscousDrag::dKidT1(const label celli) const
{
    return 0.0;
}

Foam::scalar Foam::TwoFluidFoam::viscousDrag::dKidT2(const label celli) const
{
    return 0.0;
}


void Foam::TwoFluidFoam::viscousDrag::correct()
{
    Ki_ = ((9.0/2.0)*fluid_.thermo1().mu()/(sqr(d_/2)*Cc(d_)));//*pos(fluid_.alpha() - 1e-26);
}


// ************************************************************************* //
