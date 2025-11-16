/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "condensationMonodispersion.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationMonodispersion, 0);
addToRunTimeSelectionTable(condensationModel, condensationMonodispersion, params);


condensationMonodispersion::condensationMonodispersion
(
    volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur,
    const dictionary& dict
)
:
    condensationModel(alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, satur, dict),

    n_(
        IOobject
        (
            "n",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    pNucleation_(nucleationModel::New(dict, gasThermo_, liquidThermo_, saturation_, gasProps_, liquidProps_)),
    nucleation_(pNucleation_()),

    pGrowth_(growthModel::New(dict, gasThermo_, liquidThermo_, saturation_)),
    growth_(pGrowth_())
{}

void condensationMonodispersion::constrainN()
{
    const volScalarField& J = nucleation_.J();
    forAll(alpha_, i)
    {
        if (alpha_[i] <= 1e-29 && J[i] <= 0)
        {
            n_[i] = 0;
        }
        
        n_[i] = max(n_[i], 0);
    }
}

void condensationMonodispersion::correct()
{
    // Droplet nucleation rate [m^-3 s^-1]
    /*volScalarField J
    (
        IOobject
        (
            "J",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless/dimTime/dimVolume, 0.0)
    );*/

    // Droplet growth rate [m/s]
    /*volScalarField rDot
    (
        IOobject
        (
            "rDot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    );*/

    //correctDropletRadius();

    /*volScalarField sigma = nucleation_.sigma();
    volScalarField rc = nucleation_.criticalDropletRadius(sigma);
    J = nucleation_.J(rc, sigma);
    rDot = growth_.rDot(r_);*/



    //r_ = max(pow((3.0*alpha_)/(4*pi*rho_*(n_ + dimensionedScalar("smallN", dimensionSet(-1, 0, 0, 0, 0, 0, 0), SMALL))), 1.0/3.0), rMin_);

    /*forAll(r_, i)
    {
        // mu_g - scalar eta = 1.823e-6*sqrt(T_g[i])/(1 + 673.0/T_g[i]);
        r_[i] = 0;
        Kn_[i] = 0;
        if (alpha_[i] > 1e-28 && n_[i] > 0)
        {
            r_[i] = pow((3.0*alpha_[i])/(4*pi*rho_[i]*(n_[i] + SMALL)), 1.0/3.0);
            Kn_[i] = 1.5*mu_g[i]*sqrt(Rg.value()*T_g[i])/(2*r_[i]*p[i]);
        }
    }*/

    //r_ = pos(alpha_ - 1e-28)*pow((3.0*alpha_)/(4*pi*rho_*(n_ + dimensionedScalar("smallN", dimless/dimMass, VSMALL))), 1.0/3.0);

    /*forAll(J, i)
    {
        // Dynamic viscosity and thermal conductivity
        scalar eta = 1.823e-6*sqrt(T_g[i])/(1 + 673.0/T_g[i]);
        scalar lambda_g = 7.341e-3 - 1.013e-5*T_g[i] + 1.801e-7*sqr(T_g[i]) - 9.1e-11*pow3(T_g[i]);
        //scalar lambda_g = kappa_g[i];
        
        //L[i] = p[i]*(T_g[i] - Ts[i]) + saturation_.dpsdT(Ts[i])*Ts[i]/gasProps_.rho(p[i], Ts[i]);

        // Prandtl number
        //scalar Pr = eta*Cp_g[i]/lambda_g;
        scalar Pr = 1.0;

        // Knudsen number
        scalar Kn = 0.0;
        
        // Average droplet radius
        r_[i] = 0;
        if (alpha_[i] > 1e-28 && n_[i] > 0)
        {
            r_[i] = pow((3.0*alpha_[i])/(4*pi*rho_[i]*(n_[i] + SMALL)), 1.0/3.0);
            Kn = 1.5*eta*sqrt(Rg.value()*T_g[i])/(2*r_[i]*p[i]);
        }
        Kn_[i] = Kn;

        const scalar tau = max(1.0 - T_g[i]/647.096, 0.0);
        const scalar sigma = 235.8e-3*pow(tau, 1.256)*(1 - 0.625*tau);
        //scalar sigma = liquidProps_.sigma(p[i], T_g[i]);

        scalar dH = gasProps_.Ha(p[i], T_g[i]) - gasProps_.Ha(ps[i], T_g[i]);
        scalar ds = gasProps_.S( p[i], T_g[i]) - gasProps_.S( ps[i], T_g[i]);
        scalar dG = dH - T_g[i]*ds;
        rc[i] = 2*sigma/(rhos_l[i]*dG);

        if (T_g[i] >= Ts[i])    // Evaporation
        {
            if (r_[i] <= rMin_.value())   // Complete evaporation
            {
                n_[i] = 0;
                J[i] = 0;
                rDot[i] = 0;
            }
            else                // No new droplets, only evaporation
            {
                J[i] = 0;
                rDot[i] = (lambda_g/(r_[i]*(1 + 3.18*Kn/Pr)))*((T_l[i] - T_g[i])/(rho_l[i]*L[i]));
            }
        }
        else                    // Condensation (T < Ts)
        {
            J[i] = sqrt(2*sigma/(pi*pow3(m1_.value())))*sqr(rho_g[i])/rhos_l[i]*
                exp(-beta_.value()*4*pi*sqr(rc[i])*sigma/(3*kB*T_g[i]));

            //if (r_[i] > rc[i])
            if (r_[i] > rMin_.value())
            {
                rDot[i] = (lambda_g/(r_[i]*(1 + 3.18*Kn/Pr)))*((T_l[i] - T_g[i])/(rho_l[i]*L[i]));
            }
        }
    }*/

    //correctDropletRadius();

    nucleation_.correct();
    growth_.correct(dropletRadius());
    
    {
        /*fvScalarMatrix nEqn
        (
            fvm::ddt(alpha_, rho_, n_)
            + fvm::div(alphaRhoPhi_, n_)
            ==
            J
        );*/

        fvScalarMatrix nEqn
        (
            fvm::ddt(rho_, n_)
            + fvm::div(alphaRhoPhi_, n_)
            ==
            nucleation_.J()
        );

	    nEqn.relax();
        nEqn.solve();
    }

    constrainN();
}

tmp<volScalarField> condensationMonodispersion::dropletRadius() const
{
    const scalar pi = constant::mathematical::pi;
    const dimensionedScalar alphaMin = dimensionedScalar("alphaMin", dimless, 1e-28);

    return pos(alpha_ - alphaMin)*pow((3*alpha_)/(4*pi*rho_*(n_ + dimensionedScalar("nVSMALL", dimless/dimMass, VSMALL))), 1.0/3.0);
    //return pos(alpha_ - alphaMin)*pow((3*alpha_)/(4*pi*(n_ + dimensionedScalar("nVSMALL", dimensionSet(0, -3, 0, 0, 0, 0, 0), VSMALL))), 1/3);
}

tmp<volScalarField> condensationMonodispersion::dropletDiameter() const
{
    return 2.0*max(dropletRadius(), growth_.rMin());
}

tmp<volScalarField> condensationMonodispersion::nucleationRateMassSource() const
{
    const scalar pi = constant::mathematical::pi;

    return 4.0/3.0*pi*pow3(nucleation_.rc())*nucleation_.J()*liquidThermo_.rho();
}

tmp<volScalarField> condensationMonodispersion::growthRateMassSource() const
{
    const scalar pi = constant::mathematical::pi;

    return 4*pi*rho_*n_*sqr(dropletRadius())*growth_.rDot()*liquidThermo_.rho();
    //return 4*pi*(alpha_*rho_*n_)*sqr(dropletRadius())*growth_.rDot()*liquidThermo_.rho();
}

}
}

// ************************************************************************* //
