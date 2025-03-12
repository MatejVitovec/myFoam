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

scalar Foam::slau2::massFlux
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
    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);

    const scalar aTilde   = 0.5*(aLeft + aRight);

    //const scalar sqrtUDash = sqrt(0.5*(sqr(qLeft) + sqr(qRight)));
    const scalar sqrtUDash = sqrt(0.5*(magSqr(ULeft) + magSqr(URight)));

    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar Chi = sqr(1.0 - min(1.0, (1.0/aTilde)*sqrtUDash));

    const scalar g = -max(min(MaRelLeft, 0), -1.0)*min(max(MaRelRight, 0), 1.0);

    const scalar magVnBar = (rhoLeft*mag(qLeft) + rhoRight*mag(qRight))
        /(rhoLeft + rhoRight);
        
    const scalar magVnBarPlus  = (1.0 - g)*magVnBar + g*mag(qLeft);
    const scalar magVnBarMinus = (1.0 - g)*magVnBar + g*mag(qRight);

    return 0.5*(rhoLeft*(qLeft + magVnBarPlus)
        + rhoRight*(qRight - magVnBarMinus)
        - (Chi/aTilde)*(pRight - pLeft));
}

scalar Foam::slau2::pressureFlux
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
    const scalar qLeft  = (ULeft  & normalVector);
    const scalar qRight = (URight & normalVector);

    const scalar aTilde   = 0.5*(aLeft   + aRight);
    const scalar rhoTilde = 0.5*(rhoLeft + rhoRight);

    //const scalar sqrtUDash = Foam::sqrt(0.5*(sqr(qLeft) + sqr(qRight)));
    const scalar sqrtUDash = sqrt(0.5*(magSqr(ULeft) + magSqr(URight)));

    const scalar MaRelLeft  = qLeft /aTilde;
    const scalar MaRelRight = qRight/aTilde;
       
    const scalar PPlusLeft   = ((mag(MaRelLeft)  >= 1.0) ?
        0.5*(1.0 + sign(MaRelLeft))
        : (0.25*sqr(MaRelLeft + 1.0)*(2.0 - MaRelLeft)));

    const scalar PMinusRight = ((mag(MaRelRight) >= 1.0) ?
        0.5*(1.0 - sign(MaRelRight))
        : (0.25*sqr(MaRelRight - 1.0)*(2.0 + MaRelRight)));

    return 0.5*(pLeft + pRight)
        + 0.5*(PPlusLeft - PMinusRight)*(pLeft - pRight)
        + sqrtUDash*(PPlusLeft + PMinusRight - 1.0)*rhoTilde*aTilde;
}

// ************************************************************************* //
