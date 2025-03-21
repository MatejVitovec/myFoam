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
    Foam::pVirialGas

Description
    Linearized virial equation of state expressed in terms of pressure.

    The equation of state is given by

        p = rho R T (z0 + z1 p)

    where z0 and z1 are the pseudo-virial coefficients. Note that equation corresponds to
    the standard virial equation of state 

        p = rho R' T (1 + B' p)
    
    with R' = z0 R and B' = z1/z0.

SourceFiles
    pVirialGasI.H
    pVirialGas.C

\*---------------------------------------------------------------------------*/

#ifndef pVirialGas_H
#define pVirialGas_H

#include "autoPtr.H"
#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class pVirialGas;

template<class Specie>
inline pVirialGas<Specie> operator+
(
    const pVirialGas<Specie>&,
    const pVirialGas<Specie>&
);

template<class Specie>
inline pVirialGas<Specie> operator*
(
    const scalar,
    const pVirialGas<Specie>&
);

template<class Specie>
inline pVirialGas<Specie> operator==
(
    const pVirialGas<Specie>&,
    const pVirialGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const pVirialGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class pVirialGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class pVirialGas
:
    public Specie
{
    // Private data

        //- zeroth p-virial coefficient [-]
        scalar z0_;

        //- first p-virial coefficient [1/Pa]
        scalar z1_;

public:

    // Constructors

        //- Construct from components
        inline pVirialGas
        (
            const Specie& sp,
            const scalar& z0,
            const scalar& z1
        );

        //- Construct from dictionary
        pVirialGas(const dictionary& dict);

        //- Construct as named copy
        inline pVirialGas(const word& name, const pVirialGas&);

        //- Construct and return a clone
        inline autoPtr<pVirialGas> clone() const;

        // Selector from dictionary
        inline static autoPtr<pVirialGas> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "pVirialGas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return pressure [Pa]
            inline scalar p(scalar rho, scalar T) const;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor [-]
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K)]
            inline scalar CpMCv(scalar p, scalar T) const;

// IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const pVirialGas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend pVirialGas operator+ <Specie>
        (
            const pVirialGas&,
            const pVirialGas&
        );

        friend pVirialGas operator* <Specie>
        (
            const scalar s,
            const pVirialGas&
        );

        friend pVirialGas operator== <Specie>
        (
            const pVirialGas&,
            const pVirialGas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const pVirialGas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pVirialGasI.H"

#ifdef NoRepository
    #include "pVirialGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
