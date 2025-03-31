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

#include "condensationMomentum.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"
#include "turbulentWetSteamModel.H"
#include <cassert>

#define NOHALAMA

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationMomentum, 0);
addToRunTimeSelectionTable(condensationModel, condensationMomentum, params);

condensationMomentum::condensationMomentum(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
) :
    condensationModel(rho, U, phi, gasThermo, liquidThermo),

    Q0_(
        IOobject
        (
            "Q0",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Q1_(
        IOobject
        (
            "Q1",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Q2_(
        IOobject
        (
            "Q2",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    m1_(
        "m1",
        dimMass,
        gasThermo_.lookupOrDefault<scalar>("molecularMass", 2.99046e-26)
    ),

    beta_(
        "beta",
        dimless,
        gasThermo_.lookupOrDefault<scalar>("surfaceTensionCorrection", 1.0)
    ),

    Sct_(
        "Sct",
        dimless,
        gasThermo_.lookupOrDefault<scalar>("SchmidtNumber", 0.9)
    ),

#if (OPENFOAM >= 1912)    
    kantrowitz_(
        "Kantrowitz",
        gasThermo_,
        false
    ),

    courtney_(
        "Courtney",
        gasThermo_,
        false
    )
#else
    kantrowitz_(
        gasThermo_.lookupOrDefault<bool>("Kantrowitz", false)
    ),
    courtney_(
        gasThermo_.lookupOrDefault<bool>("Courtney", false)
    )
#endif

{
    Info << "Surface tension correction beta = " << beta_.value() << nl;
    if (kantrowitz_) 
    {
        Info << "Using Kantrowitz non-isothermal correction." << nl;
    }
    if (courtney_)
    {
        Info << "Using Courtney correction." << nl;
    }

    Info << nl;
}


void condensationMomentum::correct()
{

    const scalar pi = constant::mathematical::pi;
    const scalar kB = constant::physicoChemical::k.value();
    const scalar M  = constant::physicoChemical::NA.value()*m1_.value();
    const scalar Rg = constant::physicoChemical::R.value()/M;

    const fvMesh& mesh = mesh_;
    const volScalarField w  = w();
    const volScalarField& p = liquidThermo_.p();
    const volScalarField& T = liquidThermo_.T();     //liquid temperature
    const volScalarField Ts = saturation_.Ts(p);
    const volScalarField ps = saturation_.ps(T);

    //const volScalarField sigma = steam().liquid().sigma(p, T);    //prepsano naprimo pomoci liquidProps_
    const volScalarField& rho_l = liquidThermo_.rho();
    const volScalarField& rho_g = gasThermo_.rho();
    const volScalarField L = steam().L();       // TODO latent heat implement (nejspise jako funkce ve tride condenationModel - prozatiom nejjednodusi)
    const volScalarField& Cp = gasThermo_.Cp();

    const WetSteam::turbulenceModel& turbModel =
        mesh.lookupObject<WetSteam::turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                p.group()
            )
        );
    
    const dimensionedScalar rMin("rMin", dimLength, 1e-20);
    const dimensionedScalar kg("kg", dimMass, 1.0);
    const dimensionedScalar wMin("wMin", dimless, 1.e-16);

    // Droplet nucleation rate [m^-3 s^-1]
    volScalarField J
    (
        IOobject(
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
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

    
    forAll(J, i)
    {         
        // Dynamic viscosity and thermal conductivity
        scalar eta = 1.823e-6*sqrt(T[i])/(1 + 673.0/T[i]);
        scalar lambda_g = 7.341e-3 - 1.013e-5*T[i] + 1.801e-7*sqr(T[i]) - 9.1e-11*pow3(T[i]);

        // Prandtl number
        scalar Pr = eta*Cp[i]/lambda_g;
            
        // Knudsen number
        scalar Kn = 0;

        // Average droplet radius
        scalar r = 0;
        if (w[i] > wMin.value() && Q0_[i] > 0 && Q2_[i] > sqr(rMin.value())*Q0_[i])
        {
            r = sqrt(Q2_[i]/Q0_[i]);
            Kn = 1.5*eta*sqrt(Rg*T[i])/(2*r*p[i]);
        }

        if (T[i] >= Ts[i])  // Evaporation
        {
            if (r <= rMin.value())  // Complete evaporation
            {
                w[i] = 0;
                Q0_[i] = 0;
                Q1_[i] = 0;
                Q2_[i] = 0;
                J[i] = 0;
                rDot[i] = 0;
            }
            else    // No new droplets, only evaporation
            {
                J[i] = 0;
                rDot[i] = lambda_g*(Ts[i] - T[i])/(L[i]*rho_l[i]*(1 + 3.18*Kn/Pr))/r;
            }
        }
        else    // Condensation (T < Ts)
        {
            scalar S = p[i]/ps[i];
            assert(S>=1.0);

            scalar sigma = liquidProps_.sigma(p[i], T[i]);

#ifdef HALAMA
            rc[i] = 2*sigma/(L[i]*rho_l[i]*log(Ts[i]/T[i]));
#else
            //rc[i] = 2*sigma/(rho_l[i]*Rg*T[i]*log(S));

            scalar dH = steam().gasProps().Ha(p[i],T[i]) - steam().gasProps().Ha(ps[i],T[i]);
            scalar ds = steam().gasProps().S(p[i],T[i]) - steam().gasProps().S(ps[i],T[i]);
            scalar dG = dH - T[i]*ds;
            rc[i] = 2*sigma/(rho_l[i]*dG);

#endif            
            J[i] = sqrt(2*sigma/(pi*pow3(m1_.value())))*sqr(rho_g[i])/rho_l[i]*
                exp(-beta_.value()*4*pi*sqr(rc[i])*sigma/(3*kB*T[i]));

            double gamma = Cp[i]/(Cp[i] - Rg);

            if (kantrowitz_)
            {
                scalar psi = 2*(gamma - 1)/(gamma + 1) *
                    L[i]/(Rg*T[i])*(L[i]/(Rg*T[i]) - 0.5);
                
                J[i] /= (1 + psi);      // Kantrowitz correction
            }
            if (courtney_)
            {
                J[i] /= S;              // Courtney's correction
            }

            if (r > rc[i])
            {
#ifdef HALAMA
                scalar theta = Ts[i]/T[i] - 1.0;
                scalar dTbyLog = T[i]*(1.0 + theta/2.0 - sqr(theta)/12.0 + pow3(theta)/24.0 
                    - 19.0*pow4(theta)/720.0 + 3*pow5(theta)/160.0);
                rDot[i] = lambda_g/(L[i]*rho_l[i]*(1.0 + 3.18*Kn/Pr))*
                    ((Ts[i] - T[i])/r - 2.0*sigma/(L[i]*rho_l[i]*sqr(r))*dTbyLog);
#else
                rDot[i] = lambda_g*(Ts[i] - T[i])/(L[i]*rho_l[i]*(1 + 3.18*Kn/Pr)) *
                    (r - rc[i])/sqr(r);
#endif
            }
        }
    }

    volScalarField Dt("Dt", rho_*turbModel.nut()/Sct_);
    
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
    //fields.add(w);
    fields.add(Q0_);
    fields.add(Q1_);
    fields.add(Q2_);

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi_,
            mesh.divScheme("div(phi,w_Q)")
        )
    );
    

    // Solution of Q3/w in main conservative scheme via mass source term -> alpha ~ w
    /*{
        fvScalarMatrix wEqn
        (
            fvm::ddt(rho_, w) 
            + mvConvection->fvmDiv(phi_, w) 
            - fvm::laplacian(Dt, w)
            ==
            4.0/3.0*pi*pow3(rc)*J*rho_l
            - fvm::SuSp(-4*pi*rho_*Q2_*rDot*rho_l/(w+wMin), w)
        );
        
        wEqn.relax();
        wEqn.solve();
    }*/

    {
        fvScalarMatrix Q2Eqn
        (
            fvm::ddt(rho_, Q2_) 
            + mvConvection->fvmDiv(phi_, Q2_) 
            - fvm::laplacian(Dt, Q2_)
            ==
            sqr(rc)*J + 2*rho_*Q1_*rDot
        );

        Q2Eqn.relax();
        Q2Eqn.solve();
    }

    {
        fvScalarMatrix Q1Eqn
        (
            fvm::ddt(rho_, Q1_) 
            + mvConvection->fvmDiv(phi_, Q1_) 
            - fvm::laplacian(Dt, Q1_)
            ==
            rc*J + rho_*Q0_*rDot
        );

        Q1Eqn.relax();
        Q1Eqn.solve();
    }

    {
        fvScalarMatrix Q0Eqn
        (
            fvm::ddt(rho_, Q0_) 
            + mvConvection->fvmDiv(phi_, Q0_) 
            - fvm::laplacian(Dt, Q0_)
            ==
            J
        );

        Q0Eqn.relax();
        Q0Eqn.solve();
    }
    
    //mDotGL_ = 4.0/3.0*pi*pow3(rc)*J*rho_l + 4*pi*rho_*Q2_*rDot*rho_l;

    nucleationMassSource_ = 4.0/3.0*pi*pow3(rc)*J*rho_l;
    condensationRateMassSource_ = 4*pi*rho_*Q2_*rDot*rho_l;

    // TODO
    // Explicitly constrain w & Q0 in dry steam region
    forAll(w, i)
    {
        if (w[i] <= wMin.value() && J[i] <= 0)
        {
            w[i]   = 0.0;
            Q0_[i] = 0.0;
            Q1_[i] = 0.0;
            Q2_[i] = 0.0;
        }
        w[i]   = max(w[i], 0);
        Q0_[i] = max(Q0_[i], 0);
        Q1_[i] = max(Q1_[i], 0);
        Q2_[i] = max(Q2_[i], 0);
    }


    if (mesh.time().outputTime())
    {
        J.write();
    }
}

}
}

// ************************************************************************* //
