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

#include "ausmPlusUp.H"
//#include "addToRunTimeSelectionTable.H"

/*namespace Foam
{
    defineTypeNameAndDebug(ausmPlusUp, 0);
    addToRunTimeSelectionTable(riemannSolver, ausmPlusUp, dictionary);
}*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::ausmPlusUp::ausmPlusUp(const fvMesh&, const dictionary&)
{
    /*dictionary mySubDict( dict.subOrEmptyDict("AUSMplusUpFluxCoeffs") );
    beta_  = mySubDict.lookupOrAddDefault("beta", 1.0/8.0);
    MaInf_ = mySubDict.lookupOrAddDefault("MaInf", 0.3);
    Kp_    = mySubDict.lookupOrAddDefault("Kp", 0.25);
    Ku_    = mySubDict.lookupOrAddDefault("Ku", 0.75);
    sigma_ = mySubDict.lookupOrAddDefault("sigma", 1.0);
    
    if (mySubDict.lookupOrDefault("printCoeffs", false))
        Info << mySubDict << nl;*/
}

Foam::ausmPlusUp::ausmPlusUp()
{
    beta_  = 1.0/8.0;
    MaInf_ = 0.3;
    Kp_    = 0.25;
    Ku_    = 0.75;
    sigma_ = 1.0;
}

void Foam::ausmPlusUp::calculateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rhoLeft,
    const scalar& rhoRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& pLeft,
    const scalar& pRight,
    const scalar& eLeft,
    const scalar& eRight,
    const scalar& aLeft,
    const scalar& aRight,
    const vector& Sf,
    const scalar& magSf
) const
{
    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // DensityTotalEnthalpy
    const scalar rhoHLeft  = rhoLeft *(eLeft  + 0.5*magSqr(ULeft))  + pLeft;
    const scalar rhoHRight = rhoRight*(eRight + 0.5*magSqr(URight)) + pRight;

    // DensityVelocity
    const vector rhoULeft  = rhoLeft *ULeft;
    const vector rhoURight = rhoRight*URight;

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);

    //const scalar aStar = 0.5*(aLeft + aRight);     //CRITICAL SPEED OF SOUND
    //const scalar aHatLeft  = sqr(aStar)/max(aStar,  qLeft);
    //const scalar aHatRight = sqr(aStar)/max(aStar, -qRight);
    //const scalar aTilde = min(aHatLeft, aHatRight); //COMMON SPEED OF SOUND - multiple variants possible
    const scalar aTilde = 0.5*(aLeft + aRight);
    const scalar rhoTilde = 0.5*(rhoLeft + rhoRight);
    
    const scalar sqrMaDash = (sqr(qLeft) + sqr(qRight))/(2.0*sqr(aTilde));
    const scalar sqrMaZero = min(1.0, max(sqrMaDash, sqr(MaInf_)));
    const scalar MaZero    = Foam::sqrt(sqrMaZero);

    const scalar fa = MaZero*(2.0 - MaZero);

    const scalar alpha = 3.0/16.0*(-4.0 + 5.0*sqr(fa));
    
    const scalar MaRelLeft  = qLeft/aTilde;
    const scalar MaRelRight = qRight/aTilde;
    
    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);
    
    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft  + magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight - magMaRelRight);
    
    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft  + 1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight + 1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft  - 1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight - 1.0);
    
    const scalar Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0 - 16.0*beta_*Ma2MinusLeft)));
    const scalar Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0 + 16.0*beta_*Ma2PlusRight)));
        
    const scalar Mp = -Kp_/fa*max(1.0 - sigma_*sqrMaDash, 0.0)*(pRight - pLeft)/(rhoTilde*sqr(aTilde));

    const scalar P5alphaPlusLeft   = ((magMaRelLeft >= 1.0) ?
        (Ma1PlusLeft/MaRelLeft)
        : (Ma2PlusLeft  *(( 2.0 - MaRelLeft)  - 16.0*alpha*MaRelLeft *Ma2MinusLeft)));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
        (Ma1MinusRight/MaRelRight)
        : (Ma2MinusRight*((-2.0 - MaRelRight) + 16.0*alpha*MaRelRight*Ma2PlusRight)));
    
    const scalar pU = -Ku_*P5alphaPlusLeft*P5alphaMinusRight*(rhoLeft + rhoRight)*(fa*aTilde)*(qRight-qLeft);
    
    const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight + Mp;
    const scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight + pU;
    
    const scalar URelTilde = MaRelTilde*aTilde;
    const scalar magURelTilde = mag(MaRelTilde)*aTilde;
    // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), Chap 4.2 Eq. 4.8
    // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), Chap4, Eq. 42
    rhoFlux  = (0.5*(URelTilde*(rhoLeft  + rhoRight)  - magURelTilde*(rhoRight  - rhoLeft)))*magSf;
    rhoUFlux = (0.5*(URelTilde*(rhoULeft + rhoURight) - magURelTilde*(rhoURight - rhoULeft)) + pTilde*normalVector)*magSf;
    rhoEFlux = (0.5*(URelTilde*(rhoHLeft + rhoHRight) - magURelTilde*(rhoHRight - rhoHLeft)))*magSf;
}


// ************************************************************************* //
