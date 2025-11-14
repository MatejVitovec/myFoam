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

#include "CNTNucleationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace WetSteam
{    
    defineTypeNameAndDebug(CNTNucleationModel, 0);
    addToRunTimeSelectionTable(nucleationModel, CNTNucleationModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WetSteam::CNTNucleationModel::CNTNucleationModel
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur
)
:
    nucleationModel(dict, gasThermo, liquidThermo, satur),
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

Foam::tmp<Foam::volScalarField> Foam::WetSteam::CNTNucleationModel::J(const volScalarField& rc, const volScalarField& sigma) const
{
    const volScalarField& rho_g = gasThermo_.rho();
    const volScalarField& rho_l = liquidThermo_.rho();
    const volScalarField& T_g = gasThermo_.T();
    const volScalarField& T_s = saturation_.Ts();

    const volScalarField rhos_l = saturation_.rhosl();

    const scalar pi = constant::mathematical::pi;
    const dimensionedScalar kB = constant::physicoChemical::k;

    //volScalarField J = sqrt(2*sigma/(pi*pow3(m1_)))*sqr(rho_g)/rhos_l*exp(-beta_*4*pi*sqr(rc)*sigma/(3*kB*T_g));
    //return pos(T_s - T_g)*J;

    return pos(T_s - T_g)*(sqrt(2*sigma/(pi*pow3(m1_)))*sqr(rho_g)/rhos_l*exp(-beta_*4*pi*sqr(rc)*sigma/(3*kB*T_g)));
}


// ************************************************************************* //
