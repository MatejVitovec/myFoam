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
    rhoThermo& thermo1,
    rhoThermo& thermo2
)
:
    //fluxSolver_(Foam::riemannSolver::New( 
    //    p.mesh(), 
    //    p.mesh().thisDb().lookupObject<IOdictionary>("fvSchemes")) 
    //),
    fluxSolver_(new Foam::slau2()),
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
    alpha_pos_
    (
        IOobject
        (
            "alphaL",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(alpha_))
    ),
    alpha_neg_
    (
        IOobject
        (
            "alphaR",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha_pos_
    ),
    alphaRhoFlux1_pos_
    (
        IOobject
        (
            "alphaPhi1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(alpha_*thermo1_.rho()*U1_) & mesh_.Sf())
    ),
    alphaRhoFlux1_neg_
    (
        IOobject
        (
            "alphaPhi1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(alpha_*thermo1_.rho()*U1_) & mesh_.Sf())
    ),
    alphaRhoFlux2_pos_
    (
        IOobject
        (
            "alphaPhi2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate((1.0 - alpha_)*thermo2_.rho()*U2_) & mesh_.Sf())
    ),
    alphaRhoFlux2_neg_
    (
        IOobject
        (
            "alphaPhi2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate((1.0 - alpha_)*thermo2_.rho()*U2_) & mesh_.Sf())
    ),
    alphaRhoUFlux1_pos_
    (
        IOobject
        (
            "alphaRhoUFlux1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux1_pos_*linearInterpolate(U1_)
    ),
    alphaRhoUFlux1_neg_
    (
        IOobject
        (
            "alphaRhoUFlux1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux1_neg_*linearInterpolate(U1_)
    ),
    alphaRhoUFlux2_pos_
    (
        IOobject
        (
            "alphaRhoUFlux2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux2_pos_*linearInterpolate(U2_)
    ),
    alphaRhoUFlux2_neg_
    (
        IOobject
        (
            "alphaRhoUFlux2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux2_neg_*linearInterpolate(U2_)
    ),
    alphaRhoEFlux1_pos_
    (
        IOobject
        (
            "alphaRhoEFlux1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux1_pos_*linearInterpolate(thermo1_.Cv()*T1_ + 0.5*magSqr(U1_)) // initial guess for correct dimensions
    ),
    alphaRhoEFlux1_neg_
    (
        IOobject
        (
            "alphaRhoEFlux1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux1_neg_*linearInterpolate(thermo1_.Cv()*T1_ + 0.5*magSqr(U1_)) // initial guess for correct dimensions
    ),
    alphaRhoEFlux2_pos_
    (
        IOobject
        (
            "alphaRhoEFlux2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux2_pos_*linearInterpolate(thermo2_.Cv()*T2_ + 0.5*magSqr(U2_)) // initial guess for correct dimensions
    ),
    alphaRhoEFlux2_neg_
    (
        IOobject
        (
            "alphaRhoEFlux2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRhoFlux2_neg_*linearInterpolate(thermo2_.Cv()*T2_ + 0.5*magSqr(U2_)) // initial guess for correct dimensions
    ) {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoFluidConvectiveFlux::computeFlux()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    
    surfaceScalarField pos_(IOobject("pos", mesh_), mesh_, dimensionedScalar("one", dimless, 1.0));
    surfaceScalarField neg_(IOobject("neg", mesh_), mesh_, dimensionedScalar("minusOne", dimless, -1.0));

    alpha_pos_ = interpolate(alpha_, pos_);
    alpha_neg_ = interpolate(alpha_, neg_);

    surfaceScalarField p_pos(interpolate(p_, pos_));
    surfaceScalarField p_neg(interpolate(p_, neg_));

    surfaceVectorField U1_pos(interpolate(U1_, pos_));
    surfaceVectorField U1_neg(interpolate(U1_, neg_));
    surfaceVectorField U2_pos(interpolate(U2_, pos_));
    surfaceVectorField U2_neg(interpolate(U2_, neg_));

    surfaceScalarField T1_pos(interpolate(T1_, pos_));
    surfaceScalarField T1_neg(interpolate(T1_, neg_));
    surfaceScalarField T2_pos(interpolate(T2_, pos_));
    surfaceScalarField T2_neg(interpolate(T2_, neg_));

    // Calculate fluxes at internal faces
    forAll (owner, facei)
    {
        // calculate fluxes with reconstructed primitive variables at faces
	    fluxSolver_->calculateFlux
        (
            alphaRhoFlux1_pos_[facei],  alphaRhoFlux1_neg_[facei],
            alphaRhoFlux2_pos_[facei],  alphaRhoFlux2_neg_[facei],
            alphaRhoUFlux1_pos_[facei], alphaRhoUFlux1_neg_[facei],
            alphaRhoUFlux1_pos_[facei], alphaRhoUFlux2_neg_[facei],
            alphaRhoEFlux1_pos_[facei], alphaRhoEFlux1_neg_[facei],
            alphaRhoEFlux2_pos_[facei], alphaRhoEFlux2_neg_[facei],
            alpha_pos_[facei],          alpha_neg_[facei],
            p_pos[facei],               p_neg[facei],
            U1_pos[facei],              U1_neg[facei],
            U2_pos[facei],              U2_neg[facei],
            T1_pos[facei],              T1_neg[facei],
            T2_pos[facei],              T2_neg[facei],
            Sf[facei],
            magSf[facei],
            gas1(),
            gas2()
        );
    }

    // Update boundary field and values
    forAll (alphaRhoFlux1_pos_.boundaryField(), patchi)
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

            const fvsPatchScalarField& pAlpha_pos = alpha_pos_.boundaryField()[patchi];
            const fvsPatchScalarField& pAlpha_neg = alpha_neg_.boundaryField()[patchi];
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
                fluxSolver_->calculateFlux
                (
                    pAlphaRhoFlux1_pos[facei],     pAlphaRhoFlux1_neg[facei],
                    pAlphaRhoFlux2_pos[facei],     pAlphaRhoFlux2_neg[facei],
                    pAlphaRhoUFlux1_pos[facei],    pAlphaRhoUFlux1_neg[facei],
                    pAlphaRhoUFlux2_pos[facei],    pAlphaRhoUFlux2_neg[facei],
                    pAlphaRhoEFlux1_pos[facei],    pAlphaRhoEFlux1_neg[facei],
                    pAlphaRhoEFlux2_pos[facei],    pAlphaRhoEFlux2_neg[facei],
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
                fluxSolver_->calculateFlux
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