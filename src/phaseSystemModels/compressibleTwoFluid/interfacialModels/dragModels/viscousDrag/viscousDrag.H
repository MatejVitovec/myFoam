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
    Foam::viscousDrag

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
namespace CompressibleTwoFluidFoam
{

class compressibleTwoFluid;

/*---------------------------------------------------------------------------*\
                           Class viscousDrag Declaration
\*---------------------------------------------------------------------------*/

class viscousDrag
:
    public dragModel
{
    // Private data
    
    tmp<volScalarField> Cc(const volScalarField& d) const;


public:

    //- Runtime type information
    TypeName("viscousDrag");


    // Constructors

        //- Construct from a dictionary and a compressibleTwoFluid
        viscousDrag
        (
            const dictionary& dict,
            const compressibleTwoFluid& fluid,
            const volScalarField& dropletDiameter,
            const bool registerObject
        );


    //- Destructor
    virtual ~viscousDrag();


    // Member Functions

        scalar dKidp(const label celli) const;
        scalar dKidalpha(const label celli) const;
        vector dKidU1(const label celli) const;
        vector dKidU2(const label celli) const;
        scalar dKidT1(const label celli) const;
        scalar dKidT2(const label celli) const;

        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
