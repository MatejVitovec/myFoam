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
    Foam::nucleationModel

Description

SourceFiles
    nucleationModel.C

\*---------------------------------------------------------------------------*/

#ifndef nucleationModel_H
#define nucleationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

namespace TwoFluidFoam
{

/*---------------------------------------------------------------------------*\
                           Class nucleationModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationModel
:
    public regIOobject
{
protected:

    // Protected Data
    const fluidThermo& gasThermo_;
    const fluidThermo& liquidThermo_;

    const volScalarField& sigma_;
    const volScalarField& Ts_;

    const dimensionedScalar m1_;
    const dimensionedScalar beta_;

public:

    //- Runtime type information
    TypeName("nucleationModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            nucleationModel,
            dictionary,
            (
                const dictionary& dict,
                const fluidThermo& gasThermo,
                const fluidThermo& liquidThermo,
                const volScalarField& sigma,
                const volScalarField& Ts,
                const bool registerObject
            ),
            (dict, gasThermo, liquidThermo, sigma, Ts, registerObject)
        );


    // Static data members

        //- Coefficient dimensions
        static const dimensionSet dimK;


    // Constructors

        // Construct without residual constants
        nucleationModel
        (
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const volScalarField& sigma,
            const volScalarField& Ts,
            const bool registerObject
        );

        // Construct with residual constants
        nucleationModel
        (
            const dictionary& dict,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const volScalarField& sigma,
            const volScalarField& Ts,
            const bool registerObject
        );


    //- Destructor
    virtual ~nucleationModel();


    // Selectors

        static autoPtr<nucleationModel> New
        (
            const dictionary& dict,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const volScalarField& sigma,
            const volScalarField& Ts
        );


    // Member Functions

        //- Return the rate of nucleation J
        virtual tmp<volScalarField> J(const volScalarField& rc) const;

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
