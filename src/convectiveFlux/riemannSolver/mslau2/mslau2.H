
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
    mslau2.H

\*---------------------------------------------------------------------------*/

#ifndef mslau2_H
#define mslau2_H

#include "fvCFD.H"
#include "riemannSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mslau2 Declaration
\*---------------------------------------------------------------------------*/

class mslau2
:
    public riemannSolver
{

private:
    // Private Member Functions

    //- Disallow default bitwise copy construct
    mslau2(const mslau2&) = delete;;

    //- Disallow default bitwise assignment
    void operator=(const mslau2&) = delete;

    virtual scalar massFlux
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
    ) const;

    virtual scalar pressureFlux
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
    ) const;

public:

    //TypeName("mslau2");

    // Constructors

        mslau2() : riemannSolver() {};

        mslau2(const dictionary& dict) : riemannSolver(dict) {};

    // Destructor - default
    virtual ~mslau2() {};

    // Member Functions

        // Update fluxes based on current left and right states
        void calculateFlux
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
        ) const;
        
        void calculateFlux
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
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
