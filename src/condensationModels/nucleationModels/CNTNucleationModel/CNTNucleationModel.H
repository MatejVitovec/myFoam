/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::CNTNucleationModel

Description

SourceFiles
    CNTNucleationModel.C

\*---------------------------------------------------------------------------*/

#ifndef CNTNucleationModel_H
#define CNTNucleationModel_H

#include "nucleationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                           Class CNTNucleationModel Declaration
\*---------------------------------------------------------------------------*/

class CNTNucleationModel
:
    public nucleationModel
{

private:
    // Private data

    //- Disallow default bitwise copy construct
    CNTNucleationModel(const CNTNucleationModel&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const CNTNucleationModel&) = delete;
    
protected:

    const dimensionedScalar m1_;
    const dimensionedScalar beta_;

public:

    //- Runtime type information
    TypeName("CNTNucleationModel");


    // Constructors

        // Construct from components
        CNTNucleationModel
        (
            const dictionary& dict,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const saturation& satur
        );


    //- Destructor
    virtual ~CNTNucleationModel() {}

    // Member Functions

        //- Return the rate of nucleation J
        virtual tmp<volScalarField> J(const volScalarField& rc, const volScalarField& sigma) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
