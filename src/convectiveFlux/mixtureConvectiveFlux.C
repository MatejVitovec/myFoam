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
    const volScalarField& T1,
    const volScalarField& T2,
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
    T1_(T1),
    T2_(T2),
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

    surfaceScalarField T1_pos(interpolate(T1_, pos_));
    surfaceScalarField T1_neg(interpolate(T1_, neg_));

    surfaceScalarField T2_pos(interpolate(T2_, pos_));
    surfaceScalarField T2_neg(interpolate(T2_, neg_));

    //TODO

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        const scalar rho1_pos = gasProps1().rho(p_pos[faceI], T1_pos[faceI]);
        const scalar rho1_neg = gasProps1().rho(p_neg[faceI], T1_neg[faceI]);
        const scalar rho2_pos = gasProps2().rho(p_pos[faceI], T2_pos[faceI]);
        const scalar rho2_neg = gasProps2().rho(p_neg[faceI], T2_neg[faceI]);
        
        const scalar rho_pos  = rho1_pos*(1.0 - alpha_pos[faceI]) + rho2_pos*alpha_pos[faceI];
        const scalar rho_neg  = rho1_neg*(1.0 - alpha_neg[faceI]) + rho2_neg*alpha_neg[faceI];

        const scalar w_pos = (rho2_pos/rho_pos)*alpha_pos[faceI];
        const scalar w_neg = (rho2_neg/rho_neg)*alpha_neg[faceI];

        const scalar e_pos = gasProps1().Es(p_pos[faceI], T1_pos[faceI])*(1.0 - w_pos) + gasProps2().Es(p_pos[faceI], T2_pos[faceI])*w_pos;
        const scalar e_neg = gasProps1().Es(p_neg[faceI], T1_neg[faceI])*(1.0 - w_neg) + gasProps2().Es(p_neg[faceI], T2_neg[faceI])*w_neg;

        const scalar c_pos = gasProps1().c(p_pos[faceI], T1_pos[faceI])*sqrt(1.0 - w_pos) + gasProps2().c(p_pos[faceI], T2_pos[faceI])*sqrt(w_pos);
        const scalar c_neg = gasProps1().c(p_neg[faceI], T1_neg[faceI])*sqrt(1.0 - w_neg) + gasProps2().c(p_neg[faceI], T2_neg[faceI])*sqrt(w_neg);

        // calculate fluxes with reconstructed primitive variables at faces
	    fluxSolver_->calculateFlux
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            p_pos[faceI],   p_neg[faceI],
            U_pos[faceI],   U_neg[faceI],
            T_pos[faceI],   T_neg[faceI],
            Sf[faceI],
            magSf[faceI],
            gasProps()
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
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU_pos = U_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT_pos = T_pos.boundaryField()[patchi];
            
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU_neg = U_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT_neg = T_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp_pos[facei],  pp_neg[facei],
                    pU_pos[facei],  pU_neg[facei],
                    pT_pos[facei],  pT_neg[facei],
                    pSf[facei],
                    pMagSf[facei],
                    gasProps()
                );
            }
        }
        else
        {
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU = U_.boundaryField()[patchi];
            const scalarField& pT = T_.boundaryField()[patchi];

            forAll (pp, facei)
            {
                // Calculate fluxes
                fluxSolver_->calculateFlux
                (
                    pRhoFlux[facei],
                    pRhoUFlux[facei],
                    pRhoEFlux[facei],
                    pp[facei],  pp[facei],
                    pU[facei],  pU[facei],
                    pT[facei],  pT[facei],
                    pSf[facei],
                    pMagSf[facei],
                    gasProps()
                );
            }
        }
    }
}

// ************************************************************************* //