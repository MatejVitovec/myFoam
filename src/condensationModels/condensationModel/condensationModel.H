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

    const volScalarField& alpha_;           //liquid alpha
    const volScalarField& rho_;             //liquid rho
    const volVectorField& U_;               //liquid velocity
    const surfaceScalarField& alphaRhoPhi_; //liquid alphaPhi

    const fluidThermo& gasThermo_;
    const fluidThermo& liquidThermo_;

    const saturation& saturation_;

    autoPtr<gasProperties> pGasProps_;
    const gasProperties& gasProps_;

    autoPtr<liquidProperties> pliquidProps_;
    const liquidProperties& liquidProps_;

    //autoPtr<saturationCurve> pSaturation_;
    //const saturationCurve& saturation_;

    volScalarField Kn_;

    volScalarField nucleationRateMassSource_;
    volScalarField growthRateMassSource_;
    
    // Private Member Functions

private:
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
            volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const fluidThermo& gasThermo,
            const fluidThermo& liquidThermo,
            const saturation& satur,
            const dictionary& dict
        ),
        (alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, satur, dict)
    );

    
    // Constructors

    //- Construct from components
    condensationModel
    (
        volScalarField& alpha,
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
        volScalarField& alpha,
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
    const volScalarField& nucleationRateMassSource() const
    {
        return nucleationRateMassSource_;
    }

    const volScalarField& growthRateMassSource() const
    {
        return growthRateMassSource_;
    }

    tmp<volScalarField> condensationRateMassSource() const
    {
        return nucleationRateMassSource_ + growthRateMassSource_;
    }


    const fvMesh& mesh() const
    {
        return mesh_;
    }

    /*const saturationCurve& saturation() const
    {
        return saturation_;
    }*/

    tmp<volScalarField> w() const;
    

    virtual tmp<volScalarField> dropletDiameter() const;


    //virtual tmp<volScalarField> L() const;
    
    // Check
    
    // Edit
    
    // Write
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
