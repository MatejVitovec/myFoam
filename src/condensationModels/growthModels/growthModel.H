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
    Foam::growthModel

Description

SourceFiles
    growthModel.C

\*---------------------------------------------------------------------------*/

#ifndef growthModel_H
#define growthModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluidThermo.H"
#include "saturation.H"

namespace Foam
{

namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                           Class growthModel Declaration
\*---------------------------------------------------------------------------*/

class growthModel
{
protected:

    // Protected Data

    const fluidThermo& gasThermo_;
    const fluidThermo& liquidThermo_;
    const saturation& saturation_;

    //- Minimal droplet radius
    const dimensionedScalar rMin_;

    //- Growth rate
    volScalarField rDot_;

private:

    //- Disallow default bitwise copy construct
    growthModel(const growthModel&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const growthModel&) = delete;

public:

    //- Runtime type information
    TypeName("growthModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            growthModel,
            dictionary,
            (
                const dictionary& dict,
                const fluidThermo& gasThermo,
                const fluidThermo& liquidThermo,
                const saturation& satur
            ),
            (dict, gasThermo, liquidThermo, satur)
        );


    // Constructors

        // Construct from components
        growthModel
        (
            const dictionary& dict,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const saturation& satur
        );


    //- Destructor
    virtual ~growthModel();


    // Selectors

        static autoPtr<growthModel> New
        (
            const dictionary& dict,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const saturation& satur
        );


    // Member Functions

        //- return minimal droplet radius
        dimensionedScalar rMin() const;

        //- return ther growth rate
        const volScalarField& rDot() const;

        //- correct growth rate
        virtual void correct(const volScalarField& r) = 0;

        //- correct growth rate
        virtual void correct(const volScalarField& r, const volScalarField& rc) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
