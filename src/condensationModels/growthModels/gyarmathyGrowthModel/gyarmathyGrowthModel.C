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

#include "gyarmathyGrowthModel.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace WetSteam
{    
    defineTypeNameAndDebug(gyarmathyGrowthModel, 0);
    addToRunTimeSelectionTable(growthModel, gyarmathyGrowthModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WetSteam::gyarmathyGrowthModel::gyarmathyGrowthModel
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur
)
:
    growthModel(dict, gasThermo, liquidThermo, satur),
    m1_
    (
        "m1",
        dimMass,
        dict.lookupOrDefault<scalar>("molecularMass", 2.99046e-26)
    ),
    beta_
    (
        "beta",
        dimless,
        dict.lookupOrDefault<scalar>("surfaceTensionCorrection", 1)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::WetSteam::gyarmathyGrowthModel::knudsenNumber(const volScalarField& r) const
{
    const dimensionedScalar Rg("SpecificGasConstant", dimEnergy/dimMass/dimTemperature, 461.68449); //TODO

    return 1.5*gasThermo_.mu()*sqrt(Rg*gasThermo_.T())/(2*(r + dimensionedScalar("rVSMALL", dimLength, VSMALL))*gasThermo_.p());
}

Foam::tmp<Foam::volScalarField> Foam::WetSteam::gyarmathyGrowthModel::growthRate(const volScalarField& r) const
{
    const volScalarField& T_g = gasThermo_.T();
    const volScalarField& T_l = liquidThermo_.T();
    const volScalarField& rho_l = liquidThermo_.rho();
    const volScalarField& L = saturation_.L();

    tmp<volScalarField> Kn = knudsenNumber(r);

    //const volScalarField lambda_g = gasThermo_.kappa();
    const scalar Pr = 1.0; //TODO - mozna odstranit

    //return pos(r - rMin_)*((gasThermo_.kappa()*(T_l - T_g))/(rho_l*L*(r + dimensionedScalar("rSmall", dimLength, VSMALL))*(1 + 3.18*Kn/Pr)));
    return ((gasThermo_.kappa()*(T_l - T_g))/(rho_l*L*(r + dimensionedScalar("rSmall", dimLength, VSMALL))*(1 + 3.18*Kn/Pr)));
}


void Foam::WetSteam::gyarmathyGrowthModel::correct(const volScalarField& r)
{
    rDot_ = growthRate(r);
}

// ************************************************************************* //
