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

Application
    ausmPlusUp

Description
    AUSM+up flux class, computes numerical fluxes based on AUSM+ scheme. Uses primitive 
    variables as input and gives back conservative numerical fluxes.

Author
    Matej Vitovec

SourceFiles
    ausmPlusUp.H

\*---------------------------------------------------------------------------*/

#ifndef ausmPlusUp_H
#define ausmPlusUp_H

#include "riemannSolver.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ausmPlusUp Declaration
\*---------------------------------------------------------------------------*/

class ausmPlusUp
:
    public riemannSolver
{

public:
    //TypeName("ausmPlusUp");

    // Constructors

        //- Constructor null
        ausmPlusUp(const fvMesh&, const dictionary&);
        ausmPlusUp();


    // Destructor - default
    virtual ~ausmPlusUp() {};


    // Member Functions
    
        //- Update fluxes based on current left and right states
        virtual void calculateFlux
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
        ) const;

private:

    scalar beta_;
    scalar MaInf_;
    scalar Kp_;
    scalar Ku_;
    scalar sigma_;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
