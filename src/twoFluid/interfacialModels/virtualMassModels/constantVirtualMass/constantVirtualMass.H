/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2015 OpenFOAM Foundation
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
    Foam::constantVirtualMass

Description

SourceFiles
    constantVirtualMass.C

\*---------------------------------------------------------------------------*/

#ifndef constantVirtualMass_H
#define constantVirtualMass_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "virtualMassModel.H"

namespace Foam
{

namespace TwoFluidFoam
{

/*---------------------------------------------------------------------------*\
                         Class constantVirtualMass Declaration
\*---------------------------------------------------------------------------*/

class constantVirtualMass
:
    public virtualMassModel
{
    // Private data

        //- Constant virtualMass coefficient
        const dimensionedScalar virtualMassCoeff_;

public:

    //- Runtime type information
    TypeName("constantVirtualMass");


    // Constructors

        //- Construct from a dictionary and a twoFluid
        constantVirtualMass
        (
            const twoFluid& fluid,
            const bool registerObject
        );

        //- Construct from a dictionary and a twoFluid
        constantVirtualMass
        (
            const dictionary& dict,
            const twoFluid& fluid,
            const bool registerObject
        );


    //- Destructor
    virtual ~constantVirtualMass() = default;


    // Member Functions

        //- Return the virtual mass coefficient
        virtual tmp<volScalarField> Cvm() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
