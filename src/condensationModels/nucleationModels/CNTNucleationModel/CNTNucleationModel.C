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
#include "fundamentalConstants.H"


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
    const saturation& satur,
    const gasProperties& gasProps,
    const liquidProperties& liquidProps
)
:
    nucleationModel(dict, gasThermo, liquidThermo, satur, gasProps, liquidProps),
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
    ),
    kantrowitz_
    (
	    "Kantrowitz",
	    dict,
        false
    ),
    courtney_
    (
        "Courtney",
	    dict,
	    false
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::WetSteam::CNTNucleationModel::kantrowitzCorrection() const
{
    const dimensionedScalar Rg = constant::physicoChemical::k/m1_;

    const volScalarField& gamma = gasThermo_.gamma();
    const volScalarField& T_g = gasThermo_.T();
    const volScalarField& L = saturation_.L();
    
    return 1.0/(1.0 + 2.0*(gamma - 1.0)/(gamma + 1.0)*L/(Rg*T_g)*(L/(Rg*T_g) - 0.5));    
}


Foam::tmp<Foam::volScalarField> Foam::WetSteam::CNTNucleationModel::courtneyCorrection() const
{
    const volScalarField& p = gasThermo_.p();
    const volScalarField& ps = saturation_.ps();

    return 1/(p/ps);
}

Foam::tmp<Foam::volScalarField> Foam::WetSteam::CNTNucleationModel::applyCorrections() const
{
    const fvMesh& mesh = gasThermo_.p().mesh();
    tmp<Foam::volScalarField> corrections;

    //if (kantrowitz_ || courtney_)
    //{
        corrections = tmp<Foam::volScalarField>
        (
            new Foam::volScalarField
            (
                Foam::IOobject
                (
                    "nucleationCorrections",
                    mesh.time().timeName(),
                    mesh,
                    Foam::IOobject::NO_READ,
                    Foam::IOobject::NO_WRITE
                ),
                mesh,
                Foam::dimensionedScalar("zero", Foam::dimless, 1.0)
            )
        );
    //}

    if (kantrowitz_)
    {
        corrections.ref() *= kantrowitzCorrection()();
    }

    if (courtney_)
    {
        corrections.ref() *= courtneyCorrection()();
    }

    return corrections;
}


Foam::tmp<Foam::volScalarField> Foam::WetSteam::CNTNucleationModel::nucleationRate(const volScalarField& rc, const volScalarField& sigma) const
{
    const volScalarField& rho_g = gasThermo_.rho();
    const volScalarField& T_g = gasThermo_.T();
    const volScalarField& T_s = saturation_.Ts();

    const volScalarField rhos_l = saturation_.rhosl();

    const scalar pi = constant::mathematical::pi;
    const dimensionedScalar kB = constant::physicoChemical::k;

    return pos(T_s - T_g)
        *applyCorrections()
        *(sqrt(2*sigma/(pi*pow3(m1_)))*sqr(rho_g)/rhos_l*exp(-beta_*4*pi*sqr(rc)*sigma/(3*kB*T_g)));
}

void Foam::WetSteam::CNTNucleationModel::correct()
{
    const volScalarField sigma_ = sigma();

    rc_ = criticalDropletRadius(sigma_);
    J_ = nucleationRate(rc_, sigma_);
}


// ************************************************************************* //
