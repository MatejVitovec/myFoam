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
#include "compressibleMixtureConvectiveFlux.H"
#include "directionInterpolate.H"

#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::compressibleMixtureConvectiveFlux::compressibleMixtureConvectiveFlux
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& w,
    const rhoThermo& thermo1,
    const rhoThermo& thermo2
)
:
    fluxSolver_(Foam::riemannSolver::New
    (
        p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    ),
    mesh_(p.mesh()),
    p_(p),
    U_(U),
    T_(T),
    w_(w),
    pgasProps1_(gasProperties::New(thermo1)),
    pgasProps2_(gasProperties::New(thermo2)),
    gasProps1_(pgasProps1_()),
    gasProps2_(pgasProps2_()),
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
        mesh_,
        dimensionedScalar("rhoUzero", dimDensity*dimVelocity*dimArea, 0.0)
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
        mesh_,
        dimensionedVector("rhoUUzero", dimDensity*dimVelocity*dimVelocity*dimArea, vector::zero)
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
        mesh_,
        dimensionedScalar("rhoEUzero", dimDensity*dimEnergy/dimMass*dimVelocity*dimArea, 0.0)
    ),
    rhowFlux_
    (
        IOobject
        (
            "rhowFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhowUzero", dimDensity*dimVelocity*dimArea, 0.0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::compressibleMixtureConvectiveFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    
    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    surfaceScalarField p_pos(interpolate(p_, pos_));
    surfaceScalarField p_neg(interpolate(p_, neg_));

    surfaceVectorField U_pos(interpolate(U_, pos_));
    surfaceVectorField U_neg(interpolate(U_, neg_));

    surfaceScalarField T_pos(interpolate(T_, pos_));
    surfaceScalarField T_neg(interpolate(T_, neg_));

    surfaceScalarField w_pos(interpolate(w_, pos_));
    surfaceScalarField w_neg(interpolate(w_, neg_));

    // Calculate fluxes at internal faces
    forAll (owner, facei)
    {
        const scalar rho1_pos = gasProps1().rho(p_pos[facei], T_pos[facei]);
        const scalar rho1_neg = gasProps1().rho(p_neg[facei], T_neg[facei]);
        const scalar rho2_pos = gasProps2().rho(p_pos[facei], T_pos[facei]);
        const scalar rho2_neg = gasProps2().rho(p_neg[facei], T_neg[facei]);
        

        const scalar rho_pos  = 1.0/((1.0 - w_pos[facei])/rho1_pos + w_pos[facei]/rho2_pos);
        const scalar rho_neg  = 1.0/((1.0 - w_neg[facei])/rho1_neg + w_neg[facei]/rho2_neg);

        const scalar e_pos = (1.0 - w_pos[facei])*gasProps1().Es(p_pos[facei], T_pos[facei]) + w_pos[facei]*gasProps2().Es(p_pos[facei], T_pos[facei]);
        const scalar e_neg = (1.0 - w_neg[facei])*gasProps1().Es(p_neg[facei], T_neg[facei]) + w_neg[facei]*gasProps2().Es(p_neg[facei], T_neg[facei]);

        const scalar a_pos = gasProps1_.c(p_pos[facei], T_pos[facei])*sqrt(1.0 - w_pos[facei]);
        const scalar a_neg = gasProps1_.c(p_neg[facei], T_neg[facei])*sqrt(1.0 - w_neg[facei]);

        // calculate fluxes with reconstructed primitive variables at faces
	    fluxSolver_->calculateFlux
        (
            rhoFlux_[facei],
            rhoUFlux_[facei],
            rhoEFlux_[facei],
            rhowFlux_[facei],
            p_pos[facei],   p_neg[facei],
            U_pos[facei],   U_neg[facei],
            rho_pos,        rho_neg,
            e_pos,          e_neg,
            a_pos,          a_neg,
            w_pos[facei],   w_neg[facei],
            Sf[facei],
            magSf[facei]
        );
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pRhoFlux  = rhoFlux_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pRhowFlux = rhowFlux_.boundaryFieldRef()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        if (curPatch.coupled())
        {
            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pw_pos = w_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pw_neg = w_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                const scalar prho1_pos = gasProps1().rho(pp_pos[facei], pT_pos[facei]);
                const scalar prho1_neg = gasProps1().rho(pp_neg[facei], pT_neg[facei]);
                const scalar prho2_pos = gasProps2().rho(pp_pos[facei], pT_pos[facei]);
                const scalar prho2_neg = gasProps2().rho(pp_neg[facei], pT_neg[facei]);
                
                const scalar prho_pos  = 1.0/((1.0 - pw_pos[facei])/prho1_pos + pw_pos[facei]/prho2_pos);
                const scalar prho_neg  = 1.0/((1.0 - pw_neg[facei])/prho1_neg + pw_neg[facei]/prho2_neg);

                const scalar pe_pos = (1.0 - pw_pos[facei])*gasProps1().Es(pp_pos[facei], pT_pos[facei]) + pw_pos[facei]*gasProps2().Es(pp_pos[facei], pT_pos[facei]);
                const scalar pe_neg = (1.0 - pw_neg[facei])*gasProps1().Es(pp_neg[facei], pT_neg[facei]) + pw_neg[facei]*gasProps2().Es(pp_neg[facei], pT_neg[facei]);

                const scalar pa_pos = gasProps1_.c(pp_pos[facei], pT_pos[facei])*sqrt(1.0 - pw_pos[facei]);
                const scalar pa_neg = gasProps1_.c(pp_neg[facei], pT_neg[facei])*sqrt(1.0 - pw_neg[facei]);

                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pRhowFlux[facei],
                    pp_pos[facei],  pp_neg[facei],
                    pU_pos[facei],  pU_neg[facei],
                    prho_pos,       prho_neg,
                    pe_pos,         pe_neg,
                    pa_pos,         pa_neg,
                    pw_pos[facei],  pw_neg[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
        else
        {
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU = U_.boundaryField()[patchi];
            const scalarField& pT = T_.boundaryField()[patchi];
            const scalarField& pw = w_.boundaryField()[patchi];

            forAll (pp, facei)
            {
                const scalar prho1 = gasProps1().rho(pp[facei], pT[facei]);
                const scalar prho2 = gasProps2().rho(pp[facei], pT[facei]);
                const scalar prho  = 1.0/((1.0 - pw[facei])/prho1 + pw[facei]/prho2);
                const scalar pe = (1.0 - pw[facei])*gasProps1().Es(pp[facei], pT[facei]) + pw[facei]*gasProps2().Es(pp[facei], pT[facei]);
                const scalar pa = gasProps1().c(pp[facei], pT[facei])*sqrt(1.0 - pw[facei]) + gasProps2().c(pp[facei], pT[facei])*sqrt(pw[facei]);

                // Calculate fluxes
                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pRhowFlux[facei],
                    pp[facei],  pp[facei],
                    pU[facei],  pU[facei],
                    prho,       prho,
                    pe,         pe,
                    pa,         pa,
                    pw[facei], pw[facei],
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
    }
}

// ************************************************************************* //