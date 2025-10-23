/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::H2O_IAPWS97reg1

Description
    water according to IAPWS formulation
    
Note:
    water is assumed as incompressible and therofore thermo-properties
    are calculated at saturation line only

SourceFiles
    H2O_IAPWS97reg1.C

\*---------------------------------------------------------------------------*/

#ifndef H2O_IAPWS97reg1_H
#define H2O_IAPWS97reg1_H

#include "liquidProperties.H"
#include "NSRDSfunc0.H"
#include "NSRDSfunc1.H"
#include "NSRDSfunc2.H"
#include "NSRDSfunc3.H"
#include "NSRDSfunc4.H"
#include "NSRDSfunc6.H"
#include "NSRDSfunc7.H"
#include "APIdiffCoefFunc.H"

#include "thermodynamicConstants.H"
#include "IF97.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class H2O Declaration
\*---------------------------------------------------------------------------*/

class H2O_IAPWS97reg1
:
    public liquidProperties
{
    
    class IAPWS_sigma : public thermophysicalFunction
    {
    public:
        IAPWS_sigma() {};
        IAPWS_sigma(const dictionary& dict) {};

        scalar f(scalar, scalar T) const;

        virtual void writeData(Ostream& os) const
            {};

    };

    class IAPWS_voidFunction : public thermophysicalFunction
    {
    public:
        IAPWS_voidFunction() {};
        IAPWS_voidFunction(const dictionary& dict) {};

        scalar f(scalar, scalar T) const { return 0.0; }

        virtual void writeData(Ostream& os) const
            {};

    };
    
    // Private data

        IAPWS_voidFunction rho_;
        NSRDSfunc1 pv_;
        NSRDSfunc6 hl_;
        IAPWS_voidFunction Cp_;
        IAPWS_voidFunction h_;
        IAPWS_voidFunction Cpg_;
        NSRDSfunc4 B_;
        NSRDSfunc1 mu_;
        IAPWS_voidFunction mug_;
        NSRDSfunc0 kappa_;
        IAPWS_voidFunction kappag_;
        IAPWS_sigma sigma_;
        APIdiffCoefFunc D_;


public:

    friend class liquidProperties;

    //- Runtime type information
    TypeName("H2O_IAPWS97reg1");


    // Constructors

        //- Construct null
        H2O_IAPWS97reg1();

        //- Construct from components
        H2O_IAPWS97reg1
        (
            const liquidProperties& l,
            const IAPWS_voidFunction& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc6& heatOfVapourisation,
            const IAPWS_voidFunction& heatCapacity,
            const IAPWS_voidFunction& enthalpy,
            const IAPWS_voidFunction& idealGasHeatCapacity,
            const NSRDSfunc4& secondVirialCoeff,
            const NSRDSfunc1& dynamicViscosity,
            const IAPWS_voidFunction& vapourDynamicViscosity,
            const NSRDSfunc0& thermalConductivity,
            const IAPWS_voidFunction& vapourThermalConductivity,
            const IAPWS_sigma& surfaceTension,
            const APIdiffCoefFunc& vapourDiffussivity
        );

        //- Construct from dictionary
        H2O_IAPWS97reg1(const dictionary& dict);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new H2O_IAPWS97reg1(*this));
        }


    // Member Functions

        //- Liquid density [kg/m^3]
        inline scalar rho(scalar p, scalar T) const;

        //Liquid compressibility rho/p [s^2/m^2]
        inline scalar psi(scalar p, scalar T) const;

        // Return (Cp - Cv) [J/(kg K]
        inline scalar CpMCv(scalar p, scalar T) const;

        // Entropy [J/(kg K)]
        scalar S(const scalar p, const scalar T) const;

    

        //- Vapour pressure [Pa]
        inline scalar pv(scalar p, scalar T) const;

        //- Heat of vapourisation [J/kg]
        inline scalar hl(scalar p, scalar T) const;

        //- Liquid heat capacity [J/(kg K)]
        inline scalar Cp(scalar p, scalar T) const;

        //- Liquid Enthalpy [J/(kg)]
        inline scalar h(scalar p, scalar T) const;

        //- Ideal gas heat capacity [J/(kg K)]
        inline scalar Cpg(scalar p, scalar T) const;

        //- Second Virial Coefficient [m^3/kg]
        inline scalar B(scalar p, scalar T) const;

        //- Liquid viscosity [Pa s]
        inline scalar mu(scalar p, scalar T) const;

        //- Vapour viscosity [Pa s]
        inline scalar mug(scalar p, scalar T) const;

        //- Liquid thermal conductivity  [W/(m K)]
        inline scalar kappa(scalar p, scalar T) const;

        //- Vapour thermal conductivity  [W/(m K)]
        inline scalar kappag(scalar p, scalar T) const;

        //- Surface tension [N/m]
        inline scalar sigma(scalar p, scalar T) const;

        //- Vapour diffussivity [m2/s]
        inline scalar D(scalar p, scalar T) const;

        //- Vapour diffussivity [m2/s] with specified binary pair
        inline scalar D(scalar p, scalar T, scalar Wb) const;


    // I-O

        //- Write the function coefficients
        void writeData(Ostream& os) const;

        //- Ostream Operator
        friend Ostream& operator<<(Ostream& os, const H2O_IAPWS97reg1& l);
};


Ostream& operator<<(Ostream& os, const H2O_IAPWS97reg1& l);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "H2O_IAPWS97reg1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
