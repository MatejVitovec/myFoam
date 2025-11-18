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
    Foam::functionObjects::totalLoss

Group
    grpFieldFunctionObjects

Description
    Computes the ralaxation loss as a \c volScalarField and computes the total relaxation loss.

SourceFiles
    totalLoss.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_totalLoss_H
#define functionObjects_totalLoss_H

#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class totalLoss Declaration
\*---------------------------------------------------------------------------*/

class totalLoss
:
    public fvMeshFunctionObject
{
public:

    //- Runtime type information
    TypeName("totalLoss");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        totalLoss
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        totalLoss(const totalLoss&) = delete;

        //- No copy assignment
        void operator=(const totalLoss&) = delete;


    //- Destructor
    virtual ~totalLoss() = default;


    // Member Functions

        //- Read the reaction rate data
        virtual bool read(const dictionary&);

        //- Do nothing
        virtual bool execute();

        //- Write the reaction rate fields
        virtual bool write();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

