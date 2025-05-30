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
    slau2

Description
    SLAU2 flux class, computes numerical fluxes based on SLAU2 scheme. Uses primitive 
    variables as input and gives back conservative numerical fluxes.

Author
    Matej Vitovec

SourceFiles
    slau2.H

\*---------------------------------------------------------------------------*/

#ifndef slau2_H
#define slau2_H

#include "advectionSplitting.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class slau2 Declaration
\*---------------------------------------------------------------------------*/

class slau2 : public advectionSplitting
{
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
    //TypeName("slau2");

    // Constructors
    //- Constructor null
    slau2(const fvMesh&, const dictionary&) : advectionSplitting() {};
    slau2() : advectionSplitting() {};

    // Destructor - default
    virtual ~slau2() {};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
