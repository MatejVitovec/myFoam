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
    Foam::saturation

Description

SourceFiles
    saturation.C

\*---------------------------------------------------------------------------*/

#ifndef saturation_H
#define saturation_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "fluidThermo.H"

#include "saturationCurve.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class saturation Declaration
\*---------------------------------------------------------------------------*/

class saturation
{
protected:

    // Protected Data
    const fluidThermo& thermo_;

    autoPtr<saturationCurve> pSaturation_;
    const saturationCurve& saturationCurve_;

    volScalarField ps_;
    volScalarField Ts_;

    volScalarField L_;
    
public:

    //- Runtime type information
    //TypeName("saturation");

    // Declare runtime construction    

        /*declareRunTimeSelectionTable
        (
            autoPtr,
            saturation,
            dict,
            (
                const dictionary& dict
            ),
            (dict)
        );*/

    
    // Static data members


    // Constructors

        // Construct
        saturation(const dictionary& dict, const fluidThermo& thermo);

        //saturation(const fluidThermo& thermo);

    
    
    // Selectors

        //- Select null constructed
        //static autoPtr<saturation> New(const dictionary& dict);


    //- Destructor
    virtual ~saturation() {};


    // Member Functions

    //- Saturation pressure [Pa]
    virtual const volScalarField& ps() const;
    
    //- Saturation temperature [K]
    virtual const volScalarField& Ts() const;

    //- Latent heat
    virtual const volScalarField& L() const;

    //- Derivative of saturation pressure by temperature [Pa/K]
    virtual tmp<volScalarField> dpsdT() const;

    virtual tmp<volScalarField> rhosl() const;

    virtual tmp<volScalarField> rhosv() const;

    virtual tmp<volScalarField> hsl() const;

    virtual tmp<volScalarField> hsv() const;

    virtual tmp<volScalarField> esl() const;

    virtual tmp<volScalarField> esv() const;

    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
