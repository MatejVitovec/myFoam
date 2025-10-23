/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "convectiveFlux.H"
#include "directionInterpolate.H"

#include "slau2.H"
#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::convectiveFlux::convectiveFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    psiThermo& thermo
)
:
    //fluxSolver_(Foam::riemannSolver::New( 
    //    p.mesh(), 
    //    p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    //),
    fluxSolver_(new Foam::hllc()),
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    thermo_(thermo),
    rhoFlux_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(thermo_.rho()*U_) & mesh_.Sf())
    ),
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(U_)
    ),
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(thermo_.Cv()*T_ + 0.5*magSqr(U_)) //TODO
    ) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::convectiveFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    // Thermodynamics
    const volScalarField e = thermo_.he();
    const volScalarField rho = thermo_.rho();
    
    const volScalarField a = sqrt(thermo_.Cp()/thermo_.Cv()/thermo_.psi());
    

    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    surfaceScalarField p_pos( interpolate(p_, pos_) );
    surfaceScalarField p_neg( interpolate(p_, neg_) );

    surfaceVectorField U_pos( interpolate(U_, pos_) );
    surfaceVectorField U_neg( interpolate(U_, neg_) );

    //surfaceScalarField T_pos( interpolate(T_, pos_) );
    //surfaceScalarField T_neg( interpolate(T_, neg_) );

    surfaceScalarField rho_pos = ( interpolate(rho, pos_, T_.name()) );
    surfaceScalarField rho_neg = ( interpolate(rho, neg_, T_.name()) );

    surfaceScalarField e_pos = ( interpolate(e, pos_, T_.name()) );
    surfaceScalarField e_neg = ( interpolate(e, neg_, T_.name()) );

    surfaceScalarField a_pos = ( interpolate(a, pos_, T_.name()) );
    surfaceScalarField a_neg = ( interpolate(a, neg_, T_.name()) );


    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        // calculate fluxes with reconstructed primitive variables at faces
	    fluxSolver_->calculateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            rho_pos[faceI], rho_neg[faceI],
            U_pos[faceI],   U_neg[faceI],
            p_pos[faceI],   p_neg[faceI],
            e_pos[faceI],   e_neg[faceI],
            a_pos[faceI],   a_neg[faceI],
            Sf[faceI],
            magSf[faceI]
        );
        /*fluxSolver_->calculateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            rho[own], rho[nei],
            U_[own],   U_[nei],
            p_[own],   p_[nei],
            e[own],   e[nei],
            a[own],   a[nei],
            Sf[faceI],
            magSf[faceI]
        );*/
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryFieldRef()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        if (curPatch.coupled())
        {
            /*
            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            //const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pRho_pos = rho_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pe_pos = e_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pa_pos = a_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            //const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pRho_neg = rho_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pe_neg = e_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pa_neg = a_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pRho_pos[facei],    pRho_neg[facei],
                    pU_pos[facei],      pU_neg[facei],
                    pp_pos[facei],      pp_neg[facei],
                    pe_pos[facei],      pe_neg[facei],
                    pa_pos[facei],      pa_neg[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
            */
        }
        else
        {
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU = U_.boundaryField()[patchi];
            //const scalarField& pT = T_.boundaryField()[patchi];
            const scalarField& pRho = rho.boundaryField()[patchi];
            const scalarField& pe = e.boundaryField()[patchi];
            const scalarField& pa = a.boundaryField()[patchi];

            forAll (pp, facei)
            {
                // Calculate fluxes
                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pRho[facei],    pRho[facei],
                    pU[facei],      pU[facei],
                    pp[facei],      pp[facei],
                    pe[facei],      pe[facei],
                    pa[facei],      pa[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
    }
}

// ************************************************************************* //