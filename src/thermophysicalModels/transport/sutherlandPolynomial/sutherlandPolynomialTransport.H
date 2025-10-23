/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::sutherlandPolynomialTransport

Group
    grpSpecieTransport

Description
    Transport package using Sutherland's formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu = A_s \frac{\sqrt{T}}{1 + T_s / T}
    \f]

SourceFiles
    sutherlandPolynomialTransportI.H
    sutherlandPolynomialTransport.C

\*---------------------------------------------------------------------------*/

#ifndef sutherlandPolynomialTransport_H
#define sutherlandPolynomialTransport_H

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo, int PolySize> class sutherlandPolynomialTransport;

template<class Thermo, int PolySize>
inline sutherlandPolynomialTransport<Thermo, PolySize> operator+
(
    const sutherlandPolynomialTransport<Thermo, PolySize>&,
    const sutherlandPolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
inline sutherlandPolynomialTransport<Thermo, PolySize> operator*
(
    const scalar,
    const sutherlandPolynomialTransport<Thermo, PolySize>&
);

template<class Thermo, int PolySize>
Ostream& operator<<
(
    Ostream&,
    const sutherlandPolynomialTransport<Thermo, PolySize>&
);


/*---------------------------------------------------------------------------*\
                    Class sutherlandPolynomialTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo, int PolySize=4>
class sutherlandPolynomialTransport
:
    public Thermo
{
    // Private Data

        // Sutherland's coefficients
        scalar As_, Ts_;

        //- Thermal conductivity polynomial coefficients
        Polynomial<PolySize> kappaCoeffs_;


    // Private Member Functions

        //- Coeffs name. Eg, "kappaCoeffs<10>"
        inline static word coeffsName(const char* name)
        {
            return word(name + ("Coeffs<" + std::to_string(PolySize) + '>'));
        }


        //- Calculate the Sutherland coefficients
        //  given two viscosities and temperatures
        inline void calcCoeffs
        (
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2
        );

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct from components
        inline sutherlandPolynomialTransport
        (
            const Thermo& t,
            const scalar As,
            const scalar Ts,
            const Polynomial<PolySize>& kappaCoeffs
        );

        //- Construct from two viscosities
        inline sutherlandPolynomialTransport
        (
            const Thermo& t,
            const scalar mu1, const scalar T1,
            const scalar mu2, const scalar T2,
            const Polynomial<PolySize>& kappaCoeffs
        );

        //- Construct as named copy
        inline sutherlandPolynomialTransport(const word&, const sutherlandPolynomialTransport&);

        //- Construct from dictionary
        explicit sutherlandPolynomialTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        sutherlandPolynomialTransport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<sutherlandPolynomialTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<sutherlandPolynomialTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sutherlandPolynomial<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const sutherlandPolynomialTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend sutherlandPolynomialTransport operator+ <Thermo, PolySize>
        (
            const sutherlandPolynomialTransport&,
            const sutherlandPolynomialTransport&
        );

        friend sutherlandPolynomialTransport operator* <Thermo, PolySize>
        (
            const scalar,
            const sutherlandPolynomialTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo, PolySize>
        (
            Ostream&,
            const sutherlandPolynomialTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sutherlandPolynomialTransportI.H"

#ifdef NoRepository
    #include "sutherlandPolynomialTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
