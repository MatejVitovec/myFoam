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
    advectionSplitting

Description
    Basic class for of inviscid numerical fluxes.

Author
    Matej Vitovec

SourceFiles
    advectionSplitting.C

\*---------------------------------------------------------------------------*/

#include "advectionSplitting.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::advectionSplitting::calculateFlux
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
    const vector normalVector = Sf/magSf;

    const scalar rhoLeft  = gas.rho(pLeft,  TLeft);
    const scalar rhoRight = gas.rho(pRight, TRight);

    const scalar aLeft  = gas.c(pLeft,  TLeft);
    const scalar aRight = gas.c(pRight, TRight);

    const scalar HLeft  = gas.Hs(pLeft,  TLeft)  + 0.5*magSqr(ULeft);
    const scalar HRight = gas.Hs(pRight, TRight) + 0.5*magSqr(URight);

    const scalar massFlux_ = massFlux(pLeft,    pRight,
                                      ULeft,    URight,
                                      rhoLeft,  rhoRight,
                                      aLeft,    aRight,
                                      normalVector);                               

    const scalar leftMassFlux  = 0.5*(massFlux_ + mag(massFlux_));
    const scalar rightMassFlux = 0.5*(massFlux_ - mag(massFlux_));                 

    const scalar pressureFlux_ = pressureFlux(pLeft,    pRight,
                                              ULeft,    URight,
                                              rhoLeft,  rhoRight,
                                              aLeft,    aRight,
                                              normalVector);

    rhoFlux  = (leftMassFlux       + rightMassFlux)*magSf;
    rhoUFlux = (leftMassFlux*ULeft + rightMassFlux*URight + pressureFlux_*normalVector)*magSf;
    rhoEFlux = (leftMassFlux*HLeft + rightMassFlux*HRight)*magSf;
}

void Foam::advectionSplitting::calculateFlux
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

    const scalar massFlux1 = massFlux(pLeft,    pRight,
                                      U1Left,   U1Right,
                                      rho1Left, rho1Right,
                                      aLeft,    aRight,
                                      normalVector);

    const scalar leftMassFlux1  = 0.5*(massFlux1 + mag(massFlux1));
    const scalar rightMassFlux1 = 0.5*(massFlux1 - mag(massFlux1));

    const scalar massFlux2 = massFlux(pLeft,        pRight,
                                      U2Left,       U2Right,
                                      rho2Left,     rho2Right,
                                      aLeft,        aRight,
                                      normalVector);
    
    const scalar leftMassFlux2  = 0.5*(massFlux2 + mag(massFlux2));
    const scalar rightMassFlux2 = 0.5*(massFlux2 - mag(massFlux2));

    const scalar pressureFlux1 = pressureFlux(pLeft,    pRight,
                                              U1Left,   U1Right,
                                              rho1Left, rho1Right,
                                              aLeft,    aRight,
                                              normalVector);

    const scalar pressureFlux2 = pressureFlux(pLeft,    pRight,
                                              U2Left,   U2Right,
                                              rho2Left, rho2Right,
                                              aLeft,    aRight,
                                              normalVector);

    alphaRhoFlux1Left   = (leftMassFlux1*alphaLeft + rightMassFlux1*alphaRight)*magSf;
    alphaRhoFlux1Right  = alphaRhoFlux1Left;

    const vector alphaRhoUFlux1Aux = leftMassFlux1*alphaLeft*U1Left + rightMassFlux1*alphaRight*U1Right;
    alphaRhoUFlux1Left  = (alphaRhoUFlux1Aux + pressureFlux1*alphaLeft *normalVector)*magSf;
    alphaRhoUFlux1Right = (alphaRhoUFlux1Aux + pressureFlux1*alphaRight*normalVector)*magSf;

    alphaRhoEFlux1Left  = (leftMassFlux1*alphaLeft*H1Left + rightMassFlux1*alphaRight*H1Right)*magSf;
    alphaRhoEFlux1Right = alphaRhoEFlux1Left;



    alphaRhoFlux2Left   = (leftMassFlux2*(1.0 - alphaLeft) + rightMassFlux2*(1.0 - alphaRight))*magSf;
    alphaRhoFlux2Right  = alphaRhoFlux2Left;

    const vector alphaRhoUFlux2Aux = leftMassFlux2*(1.0 - alphaLeft)*U2Left + rightMassFlux2*(1.0 - alphaRight)*U2Right;
    alphaRhoUFlux2Left  = (alphaRhoUFlux2Aux + pressureFlux2*(1.0 - alphaLeft) *normalVector)*magSf;
    alphaRhoUFlux2Right = (alphaRhoUFlux2Aux + pressureFlux2*(1.0 - alphaRight)*normalVector)*magSf;

    alphaRhoEFlux2Left  = (leftMassFlux2*(1.0 - alphaLeft)*H2Left + rightMassFlux2*(1.0 - alphaRight)*H2Right)*magSf;
    alphaRhoEFlux2Right = alphaRhoEFlux2Left;
}

// ************************************************************************* //