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

#include "slau2.H"
//#include "addToRunTimeSelectionTable.H"

/*namespace Foam
{
    defineTypeNameAndDebug(slau2, 0);
    addToRunTimeSelectionTable(riemannSolver, slau2, dictionary);
}*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::slau2::slau2(const fvMesh&, const dictionary&)
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

Foam::slau2::slau2()
{

}

void Foam::slau2::calculateFlux
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
    const vector normalVector = Sf/magSf;

    const scalar HLeft  = eLeft  + 0.5*magSqr(ULeft)  + pLeft/rhoLeft;
    const scalar HRight = eRight + 0.5*magSqr(URight) + pRight/rhoRight;

    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);

    const scalar aTilde   = 0.5*(aLeft   + aRight);
    const scalar rhoTilde = 0.5*(rhoLeft + rhoRight);

    const scalar sqrtUDash = Foam::sqrt(0.5*(sqr(qLeft) + sqr(qRight)));

    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar Chi = sqr(1.0 - min(1.0, (1.0/aTilde)*sqrtUDash));

    const scalar g = -max(min(MaRelLeft, 0), -1.0)*min(max(MaRelRight, 0), 1.0);

    const scalar magVnBar = (rhoLeft*mag(qLeft) + rhoRight*mag(qRight))
        /(rhoLeft + rhoRight);
        
    const scalar magVnBarPlus  = (1.0 - g)*magVnBar + g*mag(qLeft);
    const scalar magVnBarMinus = (1.0 - g)*magVnBar + g*mag(qRight);

    const scalar PPlusLeft   = ((mag(MaRelLeft)  >= 1.0) ?
        0.5*(1.0 + sign(MaRelLeft))
        : (0.25*sqr(MaRelLeft + 1.0)*(2.0 - MaRelLeft)));
    const scalar PMinusRight = ((mag(MaRelRight) >= 1.0) ?
        0.5*(1.0 - sign(MaRelRight))
        : (0.25*sqr(MaRelRight - 1.0)*(2.0 + MaRelRight)));

    const scalar pTilde = 0.5*(pLeft + pRight)
        + 0.5*(PPlusLeft - PMinusRight)*(pLeft - pRight)
        + sqrtUDash*(PPlusLeft + PMinusRight - 1.0)*rhoTilde*aTilde;


    const scalar mTilde = 0.5*(rhoLeft*(qLeft + magVnBarPlus)
        + rhoRight*(qRight - magVnBarMinus)
        - (Chi/aTilde)*(pRight - pLeft));

    const scalar magMTilde = mag(mTilde);

    rhoFlux  = (0.5*(mTilde + magMTilde)       + 0.5*(mTilde - magMTilde))*magSf;
    rhoUFlux = (0.5*(mTilde + magMTilde)*ULeft + 0.5*(mTilde - magMTilde)*URight + pTilde*normalVector)*magSf;
    rhoEFlux = (0.5*(mTilde + magMTilde)*HLeft + 0.5*(mTilde - magMTilde)*HRight)*magSf;
}


// ************************************************************************* //
