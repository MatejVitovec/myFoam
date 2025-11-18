/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::condensationModel

Description

SourceFiles
    condensationModel.C

\*---------------------------------------------------------------------------*/

#ifndef condensationMonodispersion_H
#define condensationMonodispersion_H

#include "condensationModel.H"
#include "nucleationModel.H"
#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------* \
                         Class condensationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationMonodispersion
:
    public condensationModel
{
private:
    // Private data

    //- Disallow default bitwise copy construct
    condensationMonodispersion(const condensationMonodispersion&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const condensationMonodispersion&) = delete;

protected:
    //- Number of droplets
    volScalarField n_;

    //- Nucelation model
    autoPtr<nucleationModel> pNucleation_;
    nucleationModel& nucleation_;

    //- Growth rate model
    autoPtr<growthModel> pGrowth_;
    growthModel& growth_;

    //- Explicitly constrain n in dry steam region
    void constrainN();

public:

    TypeName("monodispersion");
        
    // Constructors

        //- Construct from components
        condensationMonodispersion
        (
            volScalarField& alphaL,
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const saturation& satur,
            const dictionary& dict
        );
    
    //- Destructor
    virtual ~condensationMonodispersion() {};
    

    // Member Functions

        //- Main update/correction routine
        virtual void correct();

        //- Return droplet radius
        virtual tmp<volScalarField> dropletRadius() const;

        //- return droplet diameter
        virtual tmp<volScalarField> dropletDiameter() const;

        //- return condensation nucleation rate source
        virtual tmp<volScalarField> nucleationRateMassSource() const;

        //- return condensation droplet growth rate source
        virtual tmp<volScalarField> growthRateMassSource() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
