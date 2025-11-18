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

#ifndef condensationModel_H
#define condensationModel_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "gasProperties.H"
#include "liquidProperties.H"
#include "saturation.H"
#include "saturationCurve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                         Class condensationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationModel
{
protected:
    // Prprotected data

    //- Refernce to the mesh
    const fvMesh& mesh_;

    //- Time
    const Time& time_;

    const volScalarField& alphaL_;          //true volume fraction same for mixture/twoFluid/mixtureTwoFluid

    const volScalarField& alpha_;           //liquid/mixture alpha (if mixture alpha = 1)
    const volScalarField& rho_;             //liquid/mixture rho
    const volVectorField& U_;               //liquid/mixture velocity
    const surfaceScalarField& alphaRhoPhi_; //liquid/mixture alphaPhi

    const fluidThermo& gasThermo_;
    const fluidThermo& liquidThermo_;

    const saturation& saturation_;

    autoPtr<gasProperties> pGasProps_;
    const gasProperties& gasProps_;

    autoPtr<liquidProperties> pliquidProps_;
    const liquidProperties& liquidProps_;

    
private:
    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    condensationModel(const condensationModel&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const condensationModel&) = delete;


public:
    TypeName("condensationModel");

    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        condensationModel,
        params,
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
        ),
        (alphaL, alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, satur, dict)
    );

    
    // Constructors

        //- Construct from components
        condensationModel
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
       

    // Selectors
    
        //- Return a reference to the selected condensation model
        static autoPtr<condensationModel> New
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
        virtual ~condensationModel() = default;
    

    // Member Functions

        //- Main update/correction routine
        virtual void correct() = 0;
        
        // Access

        //- return condensation nucleation rate source
        virtual tmp<volScalarField> nucleationRateMassSource() const = 0;

        //- return condensation droplet growth rate source
        virtual tmp<volScalarField> growthRateMassSource() const = 0;

        //- return condensation mass rate source
        tmp<volScalarField> condensationRateMassSource() const
        {
            return nucleationRateMassSource() + growthRateMassSource();
        }

        //- return mesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- return wetness
        virtual tmp<volScalarField> w() const;

        //- Return droplet radius
        virtual tmp<volScalarField> dropletRadius() const;

        //- return droplet diameter
        virtual tmp<volScalarField> dropletDiameter() const;
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
