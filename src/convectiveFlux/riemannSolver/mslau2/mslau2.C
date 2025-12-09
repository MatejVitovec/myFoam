/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    mslau2

Description
    Basic class for of inviscid numerical fluxes.

Author
    Matej Vitovec

SourceFiles
    mslau2.C

\*---------------------------------------------------------------------------*/

#include "mslau2.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(mslau2, 0);
    addToRunTimeSelectionTable(riemannSolver, mslau2, dict);
}

Foam::mslau2::mslau2()
:
    riemannSolver(),
    epsilon_(1.0e-6) //TODO
{}

Foam::mslau2::mslau2
(
    const dictionary& dict
)
:
    riemannSolver(dict),
    epsilon_(1.0e-6) //TODO
{}

scalar Foam::mslau2::massFlux
(
    const scalar& alphaLeft,
    const scalar& alphaRight,
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

    const scalar aTilde = 0.5*(aLeft + aRight);

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

    scalar massFlux = 0.5*(rhoLeft*(qLeft + magVnBarPlus)
        + rhoRight*(qRight - magVnBarMinus));

    if(abs(alphaLeft - alphaRight) > 5*epsilon_)
    {
        return massFlux + ((max(pLeft, pRight)/min(pLeft, pRight))*(1.0 - Chi) + 1)*(pLeft - pRight)/aTilde;
    }

    return massFlux + (Chi/aTilde)*(pLeft - pRight);
}

scalar Foam::mslau2::pressureFlux
(
    const scalar& alphaLeft,
    const scalar& alphaRight,
    const scalar& pLeft,
    const scalar& pRight,
    const vector& ULeft,
    const vector& URight,
    const scalar& rhoLeft,
    const scalar& rhoRight,
    const scalar& aLeft,
    const scalar& aRight,
    const scalar& lambdaMr,
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

    scalar pressureFlux = 0.5*(pLeft + pRight)
        + 0.5*(PPlusLeft - PMinusRight)*(pLeft - pRight);

    if(abs(alphaLeft - alphaRight) > 5*epsilon_)
    {
        return pressureFlux + lambdaMr*(PPlusLeft + PMinusRight - 1.0)*rhoTilde*aTilde;
    }

    return pressureFlux + sqrtUDash*(PPlusLeft + PMinusRight - 1.0)*rhoTilde*aTilde;
}


void Foam::mslau2::calculateFlux
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
    Info << "Modified slau2 is not an option for single fluid" << endl;
}

void Foam::mslau2::calculateFlux
(
    scalar& alphaRhoFlux1Left,  scalar& alphaRhoFlux1Right,
    scalar& alphaRhoFlux2Left,  scalar& alphaRhoFlux2Right,
    vector& alphaRhoUFlux1Left, vector& alphaRhoUFlux1Right,
    vector& alphaRhoUFlux2Left, vector& alphaRhoUFlux2Right,
    scalar& alphaRhoEFlux1Left, scalar& alphaRhoEFlux1Right,
    scalar& alphaRhoEFlux2Left, scalar& alphaRhoEFlux2Right,
    const scalar& alphaLeft,    const scalar& alphaRight,
    const scalar& pLeft,        const scalar& pRight,
    const vector& U1Left,       const vector& U1Right,
    const vector& U2Left,       const vector& U2Right,
    const scalar& T1Left,       const scalar& T1Right,
    const scalar& T2Left,       const scalar& T2Right,
    const vector& Sf,
    const scalar& magSf,
    const gasProperties& gas1,
    const gasProperties& gas2
) const
{
    const vector normalVector = Sf/magSf;

    const scalar rho1Left  = gas1.rho(pLeft,  T1Left);
    const scalar rho1Right = gas1.rho(pRight, T1Right);
    const scalar rho2Left  = gas2.rho(pLeft,  T2Left);
    const scalar rho2Right = gas2.rho(pRight, T2Right);

    const scalar aLeft  = 0.5*(gas1.c(pLeft,  T1Left)  + gas2.c(pLeft,  T2Left));
    const scalar aRight = 0.5*(gas1.c(pRight, T1Right) + gas2.c(pRight, T2Right));

    //const scalar aLeft  = (alphaLeft/rho1Left   + (1.0 - alphaLeft)*rho2Left)  /(alphaLeft/(rho1Left*gas1.c(pLeft,  T1Left))     + (1.0 - alphaLeft)/(rho2Left*gas2.c(pLeft,  T2Left)));
    //const scalar aRight = (alphaRight/rho1Right + (1.0 - alphaRight)*rho2Right)/(alphaRight/(rho1Right*gas1.c(pRight,  T1Right)) + (1.0 - alphaRight)/(rho2Right*gas2.c(pRight,  T2Right)));

    const scalar H1Left  = gas1.Hs(pLeft,  T1Left)  + 0.5*magSqr(U1Left);
    const scalar H1Right = gas1.Hs(pRight, T1Right) + 0.5*magSqr(U1Right);
    const scalar H2Left  = gas2.Hs(pLeft,  T2Left)  + 0.5*magSqr(U2Left);    
    const scalar H2Right = gas2.Hs(pRight, T2Right) + 0.5*magSqr(U2Right);

    const scalar massFlux1 = massFlux(alphaLeft,    alphaRight,
                                      pLeft,        pRight,
                                      U1Left,       U1Right,
                                      rho1Left,     rho1Right,
                                      aLeft,        aRight,
                                      normalVector);

    const scalar leftMassFlux1  = 0.5*(massFlux1 + mag(massFlux1));
    const scalar rightMassFlux1 = 0.5*(massFlux1 - mag(massFlux1));

    const scalar massFlux2 = massFlux(alphaLeft,    alphaRight,
                                      pLeft,        pRight,
                                      U2Left,       U2Right,
                                      rho2Left,     rho2Right,
                                      aLeft,        aRight,
                                      normalVector);
    
    const scalar leftMassFlux2  = 0.5*(massFlux2 + mag(massFlux2));
    const scalar rightMassFlux2 = 0.5*(massFlux2 - mag(massFlux2));

    const scalar lambdaMr = abs(sqrt(0.5*(magSqr(U1Left) + magSqr(U1Right))) - sqrt(0.5*(magSqr(U2Left) + magSqr(U2Right))));

    const scalar pressureFlux1 = pressureFlux(alphaLeft,    alphaRight,
                                              pLeft,        pRight,
                                              U1Left,       U1Right,
                                              rho1Left,     rho1Right,
                                              aLeft,        aRight,
                                              lambdaMr,
                                              normalVector);

    const scalar pressureFlux2 = pressureFlux(alphaLeft,    alphaRight,
                                              pLeft,        pRight,
                                              U2Left,       U2Right,
                                              rho2Left,     rho2Right,
                                              aLeft,        aRight,
                                              lambdaMr,
                                              normalVector);


    alphaRhoFlux1Left   = (leftMassFlux1*(1.0 - alphaLeft) + rightMassFlux1*(1.0 - alphaRight))*magSf;
    alphaRhoFlux1Right  = alphaRhoFlux1Left;

    const vector alphaRhoUFlux1Aux = leftMassFlux1*(1.0 - alphaLeft)*U1Left + rightMassFlux1*(1.0 - alphaRight)*U1Right;
    alphaRhoUFlux1Left  = (alphaRhoUFlux1Aux + pressureFlux1*(1.0 - alphaLeft) *normalVector)*magSf;
    alphaRhoUFlux1Right = (alphaRhoUFlux1Aux + pressureFlux1*(1.0 - alphaRight)*normalVector)*magSf;

    alphaRhoEFlux1Left  = (leftMassFlux1*(1.0 - alphaLeft)*H1Left + rightMassFlux1*(1.0 - alphaRight)*H1Right)*magSf;
    alphaRhoEFlux1Right = alphaRhoEFlux1Left;


    alphaRhoFlux2Left   = (leftMassFlux2*alphaLeft + rightMassFlux2*alphaRight)*magSf;
    alphaRhoFlux2Right  = alphaRhoFlux2Left;

    const vector alphaRhoUFlux2Aux = leftMassFlux2*alphaLeft*U2Left + rightMassFlux2*alphaRight*U2Right;
    alphaRhoUFlux2Left  = (alphaRhoUFlux2Aux + pressureFlux2*alphaLeft *normalVector)*magSf;
    alphaRhoUFlux2Right = (alphaRhoUFlux2Aux + pressureFlux2*alphaRight*normalVector)*magSf;

    alphaRhoEFlux2Left  = (leftMassFlux2*alphaLeft*H2Left + rightMassFlux2*alphaRight*H2Right)*magSf;
    alphaRhoEFlux2Right = alphaRhoEFlux2Left;
}

// ************************************************************************* //