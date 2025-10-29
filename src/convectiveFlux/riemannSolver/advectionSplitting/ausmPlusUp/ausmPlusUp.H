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
    ausmPlusUp flux class, computes numerical fluxes based on ausmPlusUp scheme. Uses primitive 
    variables as input and gives back conservative numerical fluxes.

Author
    Matej Vitovec

SourceFiles
    ausmPlusUp.H

\*---------------------------------------------------------------------------*/

#ifndef ausmPlusUp_H
#define ausmPlusUp_H

#include "advectionSplitting.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ausmPlusUp Declaration
\*---------------------------------------------------------------------------*/

class ausmPlusUp
:
    public advectionSplitting
{

private:

    const scalar kp_;
    const scalar ku_;
    const scalar beta_;
    const scalar MaInf_;

    scalar massFlux
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

    scalar pressureFlux
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
    TypeName("ausmPlusUp");


    // Constructors

        //- Constructor null
        ausmPlusUp()
        :
            advectionSplitting(),
            kp_(0.5),
            ku_(0.5),
            beta_(1/8),    
            MaInf_(0.05)
        {};

        ausmPlusUp(const dictionary& dict)
        :
            advectionSplitting(dict),
            kp_(0.5),
            ku_(0.5),
            beta_(1/8),    
            MaInf_(0.05)
        {};

    // Destructor - default
    virtual ~ausmPlusUp() {};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
