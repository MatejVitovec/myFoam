/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

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

#include "hllc.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(hllc, 0);
    addToRunTimeSelectionTable(riemannSolver, hllc, dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllc::calculateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    bool& leftState,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
    const scalar rhoLeft,
    const scalar rhoRight,
    const scalar eLeft,
    const scalar eRight,
    const scalar aLeft,
    const scalar aRight,
    const vector Sf,
    const scalar magSf
) const
{

    // normal vector
    const vector normalVector = Sf/magSf;

    // DensityVelocity
    const vector rhoULeft  = rhoLeft *ULeft;
    const vector rhoURight = rhoRight*URight;

    // DensityTotalEnergy
    const scalar rhoELeft  = rhoLeft *(eLeft  + 0.5*magSqr(ULeft));
    const scalar rhoERight = rhoRight*(eRight + 0.5*magSqr(URight));

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);

    // Step 2: Primitive variable Riemann solver for pStar
    scalar pStar = max(0, 0.5*(pLeft + pRight)
        - 0.125*(qRight - qLeft)*(rhoLeft + rhoRight)*(aLeft + aRight));

    // Step 3: compute signal speeds for face (Lax speed estimate):
    scalar SLeft;
    scalar SRight;

    if (pStar <= pLeft)
    {
        SLeft = qLeft - aLeft;
    }        
    else
    {
        SLeft = 0.5*(qLeft + qRight) - 0.5*(aLeft + aRight);
    }
        
    if (pStar <= pRight)
    {
        SRight = qRight + aRight;
    }        
    else
    {
        SRight = 0.5*(qLeft + qRight) + 0.5*(aLeft + aRight);
    }

    const scalar SStar = (pRight - pLeft + rhoLeft*qLeft*(SLeft - qLeft) - rhoRight*qRight*(SRight - qRight))
        /stabilise(rhoLeft*(SLeft - qLeft) - rhoRight*(SRight - qRight), VSMALL);
        

    // Compute pressure in star region from the right side
    const scalar pStarRight =
        rhoRight*(qRight - SRight)*(qRight - SStar) + pRight;

    // Should be equal to the left side
    const scalar pStarLeft  =
        rhoLeft*(qLeft - SLeft)*(qLeft - SStar) + pLeft;

    // Give a warning if this is not the case
    if (mag(pStarRight - pStarLeft) > 1e-6)
    {
        Info << "mag(pStarRight - pStarLeft) > VSMALL " << endl;
    }

    // Use pStarRight for pStar, as in theory, pStarRight == pStarLeft
    //pStar = pStarRight;

    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;

    if (pos(SLeft))
    {
        // compute F_l
        convectionSpeed = qLeft;
        rhoState  = rhoLeft;
        rhoUState = rhoULeft;
        rhoEState = rhoELeft;
        pState = pLeft;
        leftState = true;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);
        //scalar omegaLeft = scalar(1.0)/(SLeft - SStar);
        pStar = pStarLeft;

        // Compute left star region
        convectionSpeed = SStar;

        rhoState  = omegaLeft* (SLeft - qLeft)*rhoLeft;
        rhoUState = omegaLeft*((SLeft - qLeft)*rhoULeft + (pStar - pLeft)*normalVector);
        rhoEState = omegaLeft*((SLeft - qLeft)*rhoELeft - pLeft*qLeft + pStar*SStar);

        pState = pStar;
        leftState = true;
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);
        //scalar omegaRight = scalar(1.0)/(SRight - SStar);
        pStar = pStarRight;

        // compute right star region
        convectionSpeed = SStar;

        rhoState  = omegaRight* (SRight - qRight)*rhoRight;
        rhoUState = omegaRight*((SRight - qRight)*rhoURight + (pStar - pRight)*normalVector);
        rhoEState = omegaRight*((SRight - qRight)*rhoERight - pRight*qRight + pStar*SStar);

        pState = pStar;
        leftState = false;
    }
    else if (neg(SRight))
    {
        // compute F_r
        convectionSpeed = qRight;
        rhoState  = rhoRight;
        rhoUState = rhoURight;
        rhoEState = rhoERight;
        pState = pRight;
        leftState = false;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    rhoUFlux = (convectionSpeed*rhoUState + pState*normalVector)*magSf;
    rhoEFlux = (convectionSpeed*(rhoEState + pState))*magSf;
}

void Foam::hllc::calculateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& TLeft,
    const scalar& TRight,
    const vector& Sf,
    const scalar& magSf,
    const gasProperties& gas
) const
{
    const scalar rhoLeft  = gas.rho(pLeft,  TLeft);
    const scalar rhoRight = gas.rho(pRight, TRight);

    const scalar eLeft  = gas.Es(pLeft,  TLeft);
    const scalar eRight = gas.Es(pRight, TRight);

    const scalar aLeft  = gas.c(pLeft,  TLeft);
    const scalar aRight = gas.c(pRight, TRight);

    bool leftState;

    calculateFlux
    (
        rhoFlux,
        rhoUFlux,
        rhoEFlux,
        leftState,
        pLeft,  pRight,
        ULeft,  URight,
        rhoLeft,rhoRight,
        eLeft,  eRight,
        aLeft,  aRight,
        Sf,
        magSf
    ); 
}


void Foam::hllc::calculateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
    const scalar rhoLeft,
    const scalar rhoRight,
    const scalar eLeft,
    const scalar eRight,
    const scalar aLeft,
    const scalar aRight,
    const vector Sf,
    const scalar magSf
) const
{
    bool leftState;

    calculateFlux
    (
        rhoFlux,
        rhoUFlux,
        rhoEFlux,
        leftState,
        pLeft,  pRight,
        ULeft,  URight,
        rhoLeft,rhoRight,
        eLeft,  eRight,
        aLeft,  aRight,
        Sf,
        magSf
    ); 
}




void Foam::hllc::calculateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    scalar& rhoQFlux,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
    const scalar rhoLeft,
    const scalar rhoRight,
    const scalar eLeft,
    const scalar eRight,
    const scalar aLeft,
    const scalar aRight,
    const scalar QLeft,
    const scalar QRight,
    const vector Sf,
    const scalar magSf
) const
{
    bool leftState;

    calculateFlux
    (
        rhoFlux,
        rhoUFlux,
        rhoEFlux,
        leftState,
        pLeft,  pRight,
        ULeft,  URight,
        rhoLeft,rhoRight,
        eLeft,  eRight,
        aLeft,  aRight,
        Sf,
        magSf
    );

    if(leftState)
    {
        rhoQFlux = rhoFlux*QLeft;
    }
    else
    {
        rhoQFlux = rhoFlux*QRight;
    }
}

// ************************************************************************* //
