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

\*---------------------------------------------------------------------------*/

#include "IAPWSIF97reg1.H"
#include "mathematicalConstants.H"
#include "../../../IF97/IF97.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::IAPWSIF97reg1<Specie>::IAPWSIF97reg1(const Specie& sp)
    :
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::IAPWSIF97reg1<Specie>::IAPWSIF97reg1
(
    const word& name,
    const IAPWSIF97reg1& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::IAPWSIF97reg1 <Specie>>
Foam::IAPWSIF97reg1<Specie>::clone() const
{
    return autoPtr<IAPWSIF97reg1<Specie>>
    (
        new IAPWSIF97reg1<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::IAPWSIF97reg1<Specie>>
Foam::IAPWSIF97reg1<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<IAPWSIF97reg1<Specie>>
    (
        new IAPWSIF97reg1<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::rho(scalar p, scalar T) const
{
    return 1/IF97::reg1::v(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::H(scalar p, scalar T) const
{
    return IF97::reg1::h(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::Cp(scalar p, scalar T) const
{
    return IF97::reg1::cp(p,T);
}

template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::E(const scalar p, const scalar T) const
{
    return IF97::reg1::u(p, T);
}

template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::Cv(scalar p, scalar T) const
{
    return IF97::reg1::cv(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::S(scalar p,scalar T) const
{
    return IF97::reg1::s(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::psi(scalar p, scalar T) const
{
    return 1.0/(p*IF97::reg1::v(p,T));
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::Z(scalar p, scalar T) const
{
    return IF97::reg1::z(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97reg1<Specie>::CpMCv(scalar p, scalar T) const
{
    return IF97::reg1::cp(p,T) - IF97::reg1::cv(p,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::IAPWSIF97reg1<Specie>::operator+=
(
    const IAPWSIF97reg1<Specie>& pg
)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::IAPWSIF97reg1<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::IAPWSIF97reg1<Specie> Foam::operator+
(
    const IAPWSIF97reg1<Specie>& pg1,
    const IAPWSIF97reg1<Specie>& pg2
)
{
    return IAPWSIF97reg1<Specie>
        (
            static_cast<const Specie&>(pg1)
            + static_cast<const Specie&>(pg2)
        );
}


template<class Specie>
Foam::IAPWSIF97reg1<Specie> Foam::operator*
(
    const scalar s,
    const IAPWSIF97reg1<Specie>& pg
)
{
    return IAPWSIF97reg1<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
Foam::IAPWSIF97reg1<Specie> Foam::operator==
(
    const IAPWSIF97reg1<Specie>& pg1,
    const IAPWSIF97reg1<Specie>& pg2
)
{
    return IAPWSIF97reg1<Specie>
        (
            static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
        );
}


// ************************************************************************* //
