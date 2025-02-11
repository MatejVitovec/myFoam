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
#include "twoFluidConvectiveFlux.H"
#include "directionInterpolate.H"

#include "slau2.H"
#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::twoFluidConvectiveFlux::twoFluidConvectiveFlux
(
    const volScalarField& p,
    const volScalarField& alpha,
    const volVectorField& U1,
    const volVectorField& U2,
    const volScalarField& T1,
    const volScalarField& T2,
    psiThermo& thermo1,
    psiThermo& thermo2
)
:
    //fluxSolver_(Foam::riemannSolver::New( 
    //    p.mesh(), 
    //    p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    //),
    fluxSolver1_(new Foam::slau2()),
    fluxSolver2_(new Foam::slau2()),
    mesh_(p.mesh()),
    p_(p),
    alpha_(alpha),
    U1_(U1),
    U2_(U2),
    T1_(T1),
    T2_(T2),
    thermo1_(thermo1),
    thermo2_(thermo2),
    gasProps1_(gasProperties::New(thermo1)),
    gasProps2_(gasProperties::New(thermo2)),
    alphaRhoFlux1_
    (
        IOobject
        (
            "alphaPhi1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(alpha1*thermo1_.rho()*U1_) & mesh_.Sf())
    ),
    alphaRhoFlux2_
    (
        IOobject
        (
            "alphaPhi2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(alpha2*thermo2_.rho()*U2_) & mesh_.Sf())
    ),
    alphaRhoUFlux1_
    (
        IOobject
        (
            "alphaRhoUFlux1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux1_*linearInterpolate(U1_)
    ),
    alphaRhoUFlux2_
    (
        IOobject
        (
            "alphaRhoUFlux2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux2_*linearInterpolate(U2_)
    ),
    alphaRhoEFlux1_
    (
        IOobject
        (
            "alphaRhoEFlux1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux1_*linearInterpolate(thermo1_.Cv()*T1_ + 0.5*magSqr(U1_)) //TODO
    ),
    alphaRhoEFlux2_
    (
        IOobject
        (
            "alphaRhoEFlux2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux2_*linearInterpolate(thermo2_.Cv()*T2_ + 0.5*magSqr(U2_)) //TODO
    ) 
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoFluidConvectiveFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    
    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    surfaceScalarField alpha_pos(interpolate(alpha_, pos_));
    surfaceScalarField alpha_neg(interpolate(alpha_, neg_));

    surfaceScalarField p_pos(interpolate(p_, pos_));
    surfaceScalarField p_neg(interpolate(p_, neg_));

    surfaceVectorField U1_pos(interpolate(U1_, pos_));
    surfaceVectorField U1_neg(interpolate(U1_, neg_));
    surfaceVectorField U2_pos(interpolate(U2_, pos_));
    surfaceVectorField U2_neg(interpolate(U2_, neg_));

    surfaceScalarField T1_pos(interpolate(T_, pos_));
    surfaceScalarField T1_neg(interpolate(T_, neg_));
    surfaceScalarField T2_pos(interpolate(T_, pos_));
    surfaceScalarField T2_neg(interpolate(T_, neg_));

    // Calculate fluxes at internal faces
    forAll (owner, faceI)
    {
        
        const label own = owner[faceI];
        const label nei = neighbour[faceI];

        // calculate fluxes with reconstructed primitive variables at faces
	    fluxSolver_->calculateFlux
        (
            alphaRhoFlux1_pos_[faceI],  alphaRhoFlux1_neg_[faceI],
            alphaRhoFlux2_pos_[faceI],  alphaRhoFlux2_neg_[faceI],
            alphaRhoUFlux1_pos_[faceI], alphaRhoUFlux1_neg_[faceI],
            alphaRhoUFlux1_pos_[faceI], alphaRhoUFlux2_neg_[faceI],
            alphaRhoEFlux1_pos_[faceI], alphaRhoEFlux1_neg_[faceI],
            alphaRhoEFlux2_pos_[faceI], alphaRhoEFlux2_neg_[faceI],
            alpha_pos[faceI],           alpha_neg[faceI],
            p_pos[faceI],               p_neg[faceI],
            U1_pos[faceI],              U1_neg[faceI],
            U2_pos[faceI],              U2_neg[faceI],
            T1_pos[faceI],              T1_neg[faceI],
            T2_pos[faceI],              T2_neg[faceI],
            Sf[faceI],
            magSf[faceI],
            gas1(),
            gas2()
        );
    }

    // Update boundary field and values
    forAll (rhoFlux_.boundaryField(), patchi)
    {
        const fvPatch& curPatch = p_.boundaryField()[patchi].patch();

        // Fluxes
        fvsPatchScalarField& pAlphaRhoFlux1_pos  = alphaRhoFlux1_pos_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaRhoFlux1_neg  = alphaRhoFlux1_neg_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pAlphaRhoUFlux1_pos = alphaRhoUFlux1_pos_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pAlphaRhoUFlux1_neg = alphaRhoUFlux1_neg_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaRhoEFlux1_pos = alphaRhoEFlux1_pos_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaRhoEFlux1_neg = alphaRhoEFlux1_neg_.boundaryFieldRef()[patchi];

        fvsPatchScalarField& pAlphaRhoFlux2_pos  = alphaRhoFlux2_pos_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaRhoFlux2_neg  = alphaRhoFlux2_neg_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pAlphaRhoUFlux2_pos = alphaRhoUFlux2_pos_.boundaryFieldRef()[patchi];
        fvsPatchVectorField& pAlphaRhoUFlux2_neg = alphaRhoUFlux2_neg_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaRhoEFlux2_pos = alphaRhoEFlux2_pos_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pAlphaRhoEFlux2_neg = alphaRhoEFlux2_neg_.boundaryFieldRef()[patchi];

        // Face areas
        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        if (curPatch.coupled())
        {
            
            // Patch fields
            const fvsPatchScalarField& pp_pos = p_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pp_neg = p_neg.boundaryField()[patchi];

            const fvsPatchScalarField& pAlpha_pos = alpha_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pAlpha_neg = alpha_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU1_pos = U1_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU1_neg = U1_neg.boundaryField()[patchi];
            const fvsPatchVectorField& pU2_pos = U2_pos.boundaryField()[patchi];
            const fvsPatchVectorField& pU2_neg = U2_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT1_pos = T1_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT1_neg = T1_neg.boundaryField()[patchi];
            const fvsPatchScalarField& pT2_pos = T2_pos.boundaryField()[patchi];
            const fvsPatchScalarField& pT2_neg = T2_neg.boundaryField()[patchi];
            
            forAll (curPatch, facei)
            {
                fluxSolver1_->calculateFlux
                (
                    pAlphaRhoFlux1_pos_[faceI],     pAlphaRhoFlux1_neg_[faceI],
                    pAlphaRhoFlux2_pos_[faceI],     pAlphaRhoFlux2_neg_[faceI],
                    pAlphaRhoUFlux1_pos_[faceI],    pAlphaRhoUFlux1_neg_[faceI],
                    pAlphaRhoUFlux2_pos_[faceI],    pAlphaRhoUFlux2_neg_[faceI],
                    pAlphaRhoEFlux1_pos_[faceI],    pAlphaRhoEFlux1_neg_[faceI],
                    pAlphaRhoEFlux2_pos_[faceI],    pAlphaRhoEFlux2_neg_[faceI],
                    pAlpha_pos[facei],              pAlpha_neg[facei],
                    pp_pos[facei],                  pp_neg[facei],
                    pU1_pos[facei],                 pU1_neg[facei],
                    pU2_pos[facei],                 pU2_neg[facei],
                    pT1_pos[facei],                 pT1_neg[facei],
                    pT2_pos[facei],                 pT2_neg[facei],
                    pSf[facei],
                    pMagSf[facei],
                    gas1(),
                    gas2()
                );
            }
        }
        else
        {
            const scalarField& pAlpha = alpha_.boundaryField()[patchi];
            const fvPatchScalarField& pp = p_.boundaryField()[patchi];
            const vectorField& pU1 = U1_.boundaryField()[patchi];
            const vectorField& pU2 = U2_.boundaryField()[patchi];
            const scalarField& pT1 = T1_.boundaryField()[patchi];
            const scalarField& pT2 = T2_.boundaryField()[patchi];

            forAll (pp, facei)
            {
                // Calculate fluxes
                fluxSolver1_->calculateFlux
                (
                    pAlphaRhoFlux1_pos[facei],      pAlphaRhoFlux1_neg[facei],
                    pAlphaRhoFlux2_pos[facei],      pAlphaRhoFlux2_neg[facei],
                    pAlphaRhoUFlux1_pos[facei],     pAlphaRhoUFlux1_neg[facei],
                    pAlphaRhoUFlux2_pos[facei],     pAlphaRhoUFlux2_neg[facei],
                    pAlphaRhoEFlux1_pos[facei],     pAlphaRhoEFlux1_neg[facei],
                    pAlphaRhoEFlux2_pos[facei],     pAlphaRhoEFlux2_neg[facei],
                    pAlpha[facei],                  pAlpha[facei],
                    pp[facei],                      pp[facei],
                    pU1[facei],                     pU1[facei],
                    pU2[facei],                     pU2[facei],
                    pT1[facei],                     pT1[facei],
                    pT2[facei],                     pT2[facei],
                    pSf[facei],
                    pMagSf[facei],
                    gas1(),
                    gas2()
                );
            }
        }
    }
}

// ************************************************************************* //