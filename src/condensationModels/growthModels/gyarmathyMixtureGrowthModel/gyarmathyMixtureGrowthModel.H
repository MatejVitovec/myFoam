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
    Foam::gyarmathyMixtureGrowthModel

Description

SourceFiles
    gyarmathyMixtureGrowthModel.C

\*---------------------------------------------------------------------------*/

#ifndef gyarmathyMixtureGrowthModel_H
#define gyarmathyMixtureGrowthModel_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                           Class gyarmathyMixtureGrowthModel Declaration
\*---------------------------------------------------------------------------*/

class gyarmathyMixtureGrowthModel
:
    public growthModel
{

private:
    // Private data

    //- Disallow default bitwise copy construct
    gyarmathyMixtureGrowthModel(const gyarmathyMixtureGrowthModel&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const gyarmathyMixtureGrowthModel&) = delete;
    
protected:

    const dimensionedScalar m1_;
    const dimensionedScalar beta_;

    //- Calculate the Knudsen number
    virtual tmp<volScalarField> knudsenNumber(const volScalarField& r) const;

    //- Calculate liquid temperature
    virtual tmp<volScalarField> liquidTemperature(const volScalarField& r, const volScalarField& rc) const;

public:

    //- Runtime type information
    TypeName("gyarmathyMixtureGrowthModel");


    // Constructors

        // Construct from components
        gyarmathyMixtureGrowthModel
        (
            const dictionary& dict,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const saturation& satur
        );


    //- Destructor
    virtual ~gyarmathyMixtureGrowthModel() {}

    // Member Functions

        //- Calculate the growth rate
        virtual tmp<volScalarField> growthRate(const volScalarField& r) const;

        //- Correct gyarmathy growth rate
        virtual void correct(const volScalarField& r);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
