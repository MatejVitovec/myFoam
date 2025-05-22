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
#include "turbulentWetSteamModel.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationMonodispersion, 0);
addToRunTimeSelectionTable(condensationModel, condensationMonodispersion, params);



condensationMonodispersion::condensationMonodispersion(
    volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const dictionary& dict
) :
    condensationModel(alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, dict),

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

    r_(
        IOobject
        (
            "r",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh()
    ),

    m1_(
        "m1",
        dimMass,
        dict.lookupOrDefault<scalar>("molecularMass", 2.99046e-26)
    ),

    beta_(
        "beta",
        dimless,
        dict.lookupOrDefault<scalar>("surfaceTensionCorrection", 1)
    ),

    Sct_(
        "Sct",
        dimless,
        dict.lookupOrDefault<scalar>("SchmidtNumber", 0.9)
    )
{
}


void condensationMonodispersion::correct()
{
    const scalar pi = constant::mathematical::pi;
    const scalar kB = constant::physicoChemical::k.value();
    const scalar M  = constant::physicoChemical::NA.value()*m1_.value();
    const scalar Rg = constant::physicoChemical::R.value()/M;

    const fvMesh& mesh = mesh_;
    const volScalarField& p = liquidThermo_.p();
    const volScalarField& T_g = gasThermo_.T();     //gas temperature
    const volScalarField& T_l = liquidThermo_.T();  //liquid temperature
    const volScalarField& rho_l = liquidThermo_.rho();
    const volScalarField& rho_g = gasThermo_.rho();

    const volScalarField Ts = saturation_.Ts(p);
    const volScalarField ps = saturation_.ps(T_g);
    const volScalarField Cp = gasThermo_.Cp();

    const scalar rMin = 1e-12;

    // Droplet nucleation rate [m^-3 s^-1]
    volScalarField J
    (
        IOobject("J", mesh.time()),
        mesh,
        dimensionedScalar("zero", dimless/dimTime/dimVolume, 0.0)
    );

    // Droplet growth rate [m/s]
    volScalarField rDot
    (
        IOobject("rDot", mesh.time()),
        mesh,
        dimensionedScalar("zero", dimLength/dimTime, 0.0)
    );

    // Droplet critical radius [m]
    volScalarField rc
    (
        IOobject("rc", mesh.time()),
        mesh,
        dimensionedScalar("zero", dimLength, 0.0)
    );

    L_ = condensationModel::L();

    forAll(J, i)
    {
        // Dynamic viscosity and thermal conductivity
        scalar eta = 1.823e-6*sqrt(T_g[i])/(1 + 673.0/T_g[i]);
        scalar lambda_g = 7.341e-3 - 1.013e-5*T_g[i] + 1.801e-7*sqr(T_g[i]) - 9.1e-11*pow3(T_g[i]);
        //L[i] = p[i]*(T_g[i] - Ts[i]) + saturation_.dpsdT(Ts[i])*Ts[i]/gasProps_.rho(p[i], Ts[i]);

        // Prandtl number
        //scalar Pr = eta*Cp[i]/lambda_g;
        scalar Pr = 0.0;

        // Knudsen number
        scalar Kn = 0;
        
        // Average droplet radius
        r[i] = 0;
        if (alpha_[i] > 1e-22 && n_[i] > 0)
        {
            r[i] = pow((3*alpha_[i])/(4*n_[i]*pi), 1/3);
            Kn = 1.5*eta*sqrt(Rg*T_g[i])/(2*r[i]*p[i]);
        }

        if (T_g[i] >= Ts[i])    // Evaporation
        {
            if (r[i] <= rMin)   // Complete evaporation
            {
                n_[i] = 0;
                J[i] = 0;
                rDot[i] = 0;
            }
            else                // No new droplets, only evaporation
            {
                J[i] = 0;
                rDot[i] = (lambda_g/(r_[i]*(1 + 3.18*Kn/Pr)))*((T_l[i] - T_g[i])/(rho_l[i]*L_[i]));
            }
        }
        else                    // Condensation (T < Ts)
        {
            scalar sigma = liquidProps_.sigma(p[i], T_g[i]);

            scalar dH = gasProps_.Ha(p[i], T_g[i]) - gasProps_.Ha(ps[i], T_g[i]);
            scalar ds = gasProps_.S( p[i], T_g[i]) - gasProps_.S( ps[i], T_g[i]);
            scalar dG = dH - T_g[i]*ds;
            rc[i] = 2*sigma/(rho_l[i]*dG);

            J[i] = sqrt(2*sigma[i]/(pi*pow3(m1_.value())))*sqr(rho_g[i])/rho_l[i]*
                exp(-beta_.value()*4*pi*sqr(rc[i])*sigma[i]/(3*kB*T_g[i]));

            if (r > rc[i])
            {
                rDot[i] = (lambda_g/(r_[i]*(1 + 3.18*Kn/Pr)))*((T_l[i] - T_g[i])/(rho_l[i]*L_[i]));
            }
        }
    }
    
    {
        fvScalarMatrix nEqn
        (
            fvm::ddt(alpha_, rho_, n_)
            + fvm::div(alphaRhoPhi_, n_)
            ==
            alpha_*rho_*J
        );

	    nEqn.relax();
        nEqn.solve();
    }

    r_ = pow((3*alpha_)/(4*n_*pi), 1/3);

    nucleationRateMassSource_ = 4.0/3.0*pi*pow3(rc)*J*rho_l;
    growthRateMassSource_ = 4*pi*(alpha_*rho_*n_)*sqr(r_)*rDot*rho_l;

    // Explicitly constrain n in dry steam region
    forAll(alpha_, i)
    {
        if (alpha_[i] <= 1e-24 && J[i] <= 0)
        {
            n_[i] = 0;
        }
        
        n_[i] = max(n_[i], 0);
    }
}

tmp<volScalarField> condensationMonodispersion::dropletDiameter() const
{
    return 2*max(r_, rMin_);
}

tmp<volScalarField> condensationMonodispersion::L() const
{
    return L_;
}


}
}

// ************************************************************************* //
