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
    Foam::dragModels::viscousDrag

Description
    Schiller and Naumann drag model for dispersed bubbly flows.

SourceFiles
    viscousDrag.C

\*---------------------------------------------------------------------------*/

#ifndef viscousDrag_H
#define viscousDrag_H

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
                           Class viscousDrag Declaration
\*---------------------------------------------------------------------------*/

class viscousDrag
:
    public dragModel
{
    // Private data

        //- Constant drag coefficient
        const dimensionedScalar dragCoeff_;

    
    tmp<volScalarField> Cc() const;


public:

    //- Runtime type information
    TypeName("viscousDrag");


    // Constructors

        //- Construct from a dictionary and a twoFluid
        viscousDrag
        (
            const dictionary& dict,
            const twoFluid& fluid,
            const bool registerObject
        );


    //- Destructor
    virtual ~viscousDrag();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> CdRe() const;

        virtual tmp<volScalarField> Ki(const volScalarField& d) const;

        virtual tmp<volScalarField> K(const volScalarField& d) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
