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
#include "mixtureConvectiveFlux.H"
#include "directionInterpolate.H"

#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::mixtureConvectiveFlux::mixtureConvectiveFlux
(
    const volScalarField& alpha,
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const psiThermo& thermo
)
:
    fluxSolver_(Foam::riemannSolver::New
    (
        p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    ),
    mesh_(p.mesh()),
    alpha_(alpha),
    p_(p),
    U_(U),
    T_(T),
    gasProps1_(gasProperties::New(thermo1)),
    gasProps2_(gasProperties::New(thermo2)),
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
    ) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixtureConvectiveFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    
    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    surfaceScalarField alpha_pos(interpolate(alpha_, pos_));
    surfaceScalarField alpha_neg(interpolate(alpha_, neg_));

    surfaceScalarField p_pos(interpolate(p_, pos_));
    surfaceScalarField p_neg(interpolate(p_, neg_));

    surfaceVectorField U_pos(interpolate(U_, pos_));
    surfaceVectorField U_neg(interpolate(U_, neg_));

    surfaceScalarField T_pos(interpolate(T_, pos_));
    surfaceScalarField T_neg(interpolate(T_, neg_));

    //TODO

    // Calculate fluxes at internal faces
    forAll (owner, facei)
    {
        const scalar rho1_pos = gasProps1().rho(p_pos[facei], T_pos[facei]);
        const scalar rho1_neg = gasProps1().rho(p_neg[facei], T_neg[facei]);
        const scalar rho2_pos = gasProps2().rho(p_pos[facei], T_pos[facei]);
        const scalar rho2_neg = gasProps2().rho(p_neg[facei], T_neg[facei]);
        
        const scalar rho_pos  = rho1_pos*(1.0 - alpha_pos[facei]) + rho2_pos*alpha_pos[facei];
        const scalar rho_neg  = rho1_neg*(1.0 - alpha_neg[facei]) + rho2_neg*alpha_neg[facei];

        const scalar w_pos = (rho2_pos/rho_pos)*alpha_pos[facei];
        const scalar w_neg = (rho2_neg/rho_neg)*alpha_neg[facei];

        const scalar e_pos = gasProps1().Es(p_pos[facei], T_pos[facei])*(1.0 - w_pos) + gasProps2().Es(p_pos[facei], T_pos[facei])*w_pos;
        const scalar e_neg = gasProps1().Es(p_neg[facei], T_neg[facei])*(1.0 - w_neg) + gasProps2().Es(p_neg[facei], T_neg[facei])*w_neg;

        const scalar a_pos = gasProps1().c(p_pos[facei], T_pos[facei])*sqrt(1.0 - w_pos) + gasProps2().c(p_pos[facei], T_pos[facei])*sqrt(w_pos);
        const scalar a_neg = gasProps1().c(p_neg[facei], T_neg[facei])*sqrt(1.0 - w_neg) + gasProps2().c(p_neg[facei], T_neg[facei])*sqrt(w_neg);

        // calculate fluxes with reconstructed primitive variables at faces
	    fluxSolver_->calculateFlux
        (
            rhoFlux_[facei],
            rhoUFlux_[facei],
            rhoEFlux_[facei],
            p_pos[facei],   p_neg[facei],
            U_pos[facei],   U_neg[facei],
            rho_pos,        rho_neg,
            e_pos,          e_neg,
            a_pos,          a_neg,
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

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        if (curPatch.coupled())
        {
            // Patch fields
            const fvsPatchScalarField& palpha_pos = alpha_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& palpha_neg = alpha_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                const scalar rho1_pos = gasProps1().rho(pp_pos[facei], pT_pos[facei]);
                const scalar rho1_neg = gasProps1().rho(pp_neg[facei], pT_neg[facei]);
                const scalar rho2_pos = gasProps2().rho(pp_pos[facei], pT_pos[facei]);
                const scalar rho2_neg = gasProps2().rho(pp_neg[facei], pT_neg[facei]);
                
                const scalar rho_pos  = rho1_pos*(1.0 - palpha_pos[facei]) + rho2_pos*palpha_pos[facei];
                const scalar rho_neg  = rho1_neg*(1.0 - palpha_neg[facei]) + rho2_neg*palpha_neg[facei];

                const scalar w_pos = (rho2_pos/rho_pos)*palpha_pos[facei];
                const scalar w_neg = (rho2_neg/rho_neg)*palpha_neg[facei];

                const scalar e_pos = gasProps1().Es(pp_pos[facei], pT_pos[facei])*(1.0 - w_pos) + gasProps2().Es(pp_pos[facei], pT_pos[facei])*w_pos;
                const scalar e_neg = gasProps1().Es(pp_neg[facei], pT_neg[facei])*(1.0 - w_neg) + gasProps2().Es(pp_neg[facei], pT_neg[facei])*w_neg;

                const scalar a_pos = gasProps1().c(pp_pos[facei], pT_pos[facei])*sqrt(1.0 - w_pos) + gasProps2().c(pp_pos[facei], pT_pos[facei])*sqrt(w_pos);
                const scalar a_neg = gasProps1().c(pp_neg[facei], pT_neg[facei])*sqrt(1.0 - w_neg) + gasProps2().c(pp_neg[facei], pT_neg[facei])*sqrt(w_neg);

                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp_pos[facei],  pp_neg[facei],
                    pU_pos[facei],  pU_neg[facei],
                    rho_pos,        rho_neg,
                    e_pos,          e_neg,
                    a_pos,          a_neg,
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
        else
        {
            const scalarField& palpha = alpha_.boundaryField()[patchi];
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU = U_.boundaryField()[patchi];
            const scalarField& pT = T_.boundaryField()[patchi];

            forAll (pp, facei)
            {
                const scalar prho1 = gasProps1().rho(pp[facei], pT[facei]);
                const scalar prho2 = gasProps2().rho(pp[facei], pT[facei]);
                const scalar prho  = prho1*(1.0 - palpha[facei]) + prho2*palpha[facei];
                const scalar pw = (prho2/prho)*palpha[facei];
                const scalar pe = gasProps1().Es(pp[facei], pT[facei])*(1.0 - pw) + gasProps2().Es(pp[facei], pT[facei])*pw;
                const scalar pa = gasProps1().c(pp[facei], pT[facei])*sqrt(1.0 - pw) + gasProps2().c(pp[facei], pT[facei])*sqrt(pw);

                // Calculate fluxes
                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp[facei],  pp[facei],
                    pU[facei],  pU[facei],
                    prho,       prho,
                    pe,         pe,
                    pa,         pa,
                    pSf[facei],
                    pMagSf[facei]
                );
            }
        }
    }
}

// ************************************************************************* //