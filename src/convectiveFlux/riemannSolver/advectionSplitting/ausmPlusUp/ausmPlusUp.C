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
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(ausmPlusUp, 0);
    addToRunTimeSelectionTable(riemannSolver, ausmPlusUp, dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar Foam::ausmPlusUp::massFlux
(
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& rhoLeft,
    const scalar& rhoRight,
    const scalar& aLeft,
    const scalar& aRight,
    const vector& normalVector
) const
{
    const scalar aTilde   = 0.5*(aLeft   + aRight);
    const scalar rhoTilde = 0.5*(rhoLeft + rhoRight);

    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);
    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;  

    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);

    //const scalar M4Plus  = (magMaRelLeft >= 1.0)  ? 0.5*(MaRelLeft  + magMaRelLeft)  : 0.25*sqr(MaRelLeft  + 1.0)*(1.0 - 16.0*beta_*0.25*sqr(MaRelLeft  - 1.0));
    //const scalar M4Minus = (magMaRelRight >= 1.0) ? 0.5*(MaRelRight - magMaRelRight) : 0.25*sqr(MaRelRight - 1.0)*(1.0 + 16.0*beta_*0.25*sqr(MaRelRight + 1.0));
    const scalar M4Plus  = (magMaRelLeft >= 1.0)  ? 0.5*(MaRelLeft  + magMaRelLeft)  :  0.25*sqr(MaRelLeft  + 1.0) - beta_*sqr(sqr(MaRelLeft) - 1.0);
    const scalar M4Minus = (magMaRelRight >= 1.0) ? 0.5*(MaRelRight - magMaRelRight) : -0.25*sqr(MaRelRight - 1.0) - beta_*sqr(sqr(MaRelRight) - 1.0);

    //const scalar MaBarSqr = (magSqr(ULeft) + magSqr(URight))/(2.0*sqr(aTilde));
    //const scalar MaBarSqr = (sqr(qLeft) + sqr(qRight))/(2.0*sqr(aTilde));
    const scalar MaBarSqr = 0.5*(sqr(MaRelLeft) + sqr(MaRelRight));
    //const scalar MaZero = sqrt(min(1.0, max(MaBarSqr, MaInf_)));
    //const scalar fa = MaZero*(2.0 - MaZero);
    const scalar MaP = kp_*(M4Plus - M4Minus + 0.5*(MaRelRight - MaRelLeft - magMaRelLeft - magMaRelRight))*max(1.0 - MaBarSqr, 0.0)*((pLeft - pRight)/(rhoTilde*sqr(aTilde)));

    const scalar Ma = M4Plus + M4Minus + MaP;

    return (Ma > 0) ? Ma*aTilde*rhoLeft : Ma*aTilde*rhoRight;
}

scalar Foam::ausmPlusUp::pressureFlux
(
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& rhoLeft,
    const scalar& rhoRight,
    const scalar& aLeft,
    const scalar& aRight,
    const vector& normalVector
) const
{
    const scalar aTilde   = 0.5*(aLeft   + aRight);
    const scalar rhoTilde = 0.5*(rhoLeft + rhoRight);

    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);
    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;  

    /*const scalar MaBarSqr = (sqr(qLeft) + sqr(qRight))/(2.0*sqr(aTilde));
    const scalar MaZero = sqrt(min(1.0, max(MaBarSqr, MaInf_)));
    const scalar fa = MaZero*(2.0 - MaZero);
    const scalar theta = (3.0/16.0)*(-4.0 + 5.0*sqr(fa));
    const scalar M2Plus  = 0.25*sqr(MaRelLeft  + 1.0);
    const scalar M2Minus = 0.25*sqr(MaRelRight - 1.0);
    const scalar PPlusLeft   = (mag(MaRelLeft)  >= 1.0) ? 0.5*(1.0 + sign(MaRelLeft))  : (1.0/MaRelLeft)*M2Plus  *(( 2.0 - MaRelLeft)  - 16.0*theta*MaRelLeft*M2Plus);
    const scalar PMinusRight = (mag(MaRelRight) >= 1.0) ? 0.5*(1.0 - sign(MaRelRight)) : (1.0/MaRelRight)*M2Minus*((-2.0 - MaRelRight) + 16.0*theta*MaRelRight*M2Minus);*/
    const scalar PPlusLeft   = (mag(MaRelLeft)  >= 1.0) ? 0.5*(1.0 + sign(MaRelLeft))  : 0.25*sqr(MaRelLeft + 1.0)*(2.0 - MaRelLeft) + (3/16)*MaRelLeft*sqr(sqr(MaRelLeft - 1.0));
    const scalar PMinusRight = (mag(MaRelRight) >= 1.0) ? 0.5*(1.0 - sign(MaRelRight)) : 0.25*sqr(MaRelRight - 1.0)*(2.0 + MaRelRight) - (3/16)*MaRelRight*sqr(sqr(MaRelRight - 1.0));

    const scalar Pu = ku_*PPlusLeft*PMinusRight*(rhoTilde*aTilde)*(qLeft - qRight);

    return PPlusLeft*pLeft + PMinusRight*pRight + Pu;
}

// ************************************************************************* //
