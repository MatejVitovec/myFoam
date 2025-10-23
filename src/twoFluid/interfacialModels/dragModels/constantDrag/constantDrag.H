/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::dragModels::constantDrag

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    constantDrag.C

\*---------------------------------------------------------------------------*/

#ifndef constantDrag_H
#define constantDrag_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{

class twoFluid;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class constantDrag Declaration
\*---------------------------------------------------------------------------*/

class constantDrag
:
    public dragModel
{
    // Private data

        //- Constant drag coefficient
        const dimensionedScalar dragCoeff_;


public:

    //- Runtime type information
    TypeName("constantDrag");


    // Constructors

        //- Construct from a dictionary and a twoFluid
        constantDrag
        (
            const dictionary& dict,
            const twoFluid& fluid,
            const bool registerObject
        );


    //- Destructor
    virtual ~constantDrag();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
