/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::totalPressure

Group
    grpFieldFunctionObjects

Description
    Computes the total pressure as a \c volScalarField and computes the total pressure.

SourceFiles
    totalPressure.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_totalPressure_H
#define functionObjects_totalPressure_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"
#include "dictionary.H"
#include "wordList.H"

#include "fluidThermo.H"
#include "gasProperties.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class totalPressure Declaration
\*---------------------------------------------------------------------------*/

class totalPressure
:
    public fvMeshFunctionObject
{
    // Field names
    word phaseName_;
    word pName_;
    word TName_;
    word UName_;
    word rhoName_;

    // Patches for averaging
    wordList patches_;

    // Write local p0 field
    bool writeField_;

    // Calculate patch average
    /*void calculatePatchAverage
    (
        const volScalarField& p0
    ) const;*/

public:

    //- Runtime type information
    TypeName("totalPressure");

    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        totalPressure
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    // Destructor
    virtual ~totalPressure() = default;

    // Member Functions

        virtual bool read(const dictionary& dict);

        virtual bool execute();

        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

