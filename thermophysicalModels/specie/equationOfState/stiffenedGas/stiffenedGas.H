/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::stiffenedGas

Description
    Stiffened gas equation of state

SourceFiles
    stiffenedGasI.H
    stiffenedGas.C

\*---------------------------------------------------------------------------*/

#ifndef stiffenedGas_H
#define stiffenedGas_H

#include "autoPtr.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class stiffenedGas;

template<class Specie>
inline stiffenedGas<Specie> operator+
(
    const stiffenedGas<Specie>&,
    const stiffenedGas<Specie>&
);

template<class Specie>
inline stiffenedGas<Specie> operator*
(
    const scalar,
    const stiffenedGas<Specie>&
);

template<class Specie>
inline stiffenedGas<Specie> operator==
(
    const stiffenedGas<Specie>&,
    const stiffenedGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const stiffenedGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class stiffenedGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class stiffenedGas
:
    public Specie
{
    // Private data

        //- pressure shift
        scalar pInf_;


public:

    // Constructors

        //- Construct from components
        inline stiffenedGas
        (
            const Specie& sp,
            const scalar& pInf
        );

        //- Construct from dictionary
        stiffenedGas(const dictionary& dict);

        //- Construct as named copy
        inline stiffenedGas(const word& name, const stiffenedGas&);

        //- Construct and return a clone
        inline autoPtr<stiffenedGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<stiffenedGas> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "stiffenedGas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(scalar p, scalar T) const;


// IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const stiffenedGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend stiffenedGas operator+ <Specie>
        (
            const stiffenedGas&,
            const stiffenedGas&
        );

        friend stiffenedGas operator* <Specie>
        (
            const scalar s,
            const stiffenedGas&
        );

        friend stiffenedGas operator== <Specie>
        (
            const stiffenedGas&,
            const stiffenedGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const stiffenedGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stiffenedGasI.H"

#ifdef NoRepository
    #include "stiffenedGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
