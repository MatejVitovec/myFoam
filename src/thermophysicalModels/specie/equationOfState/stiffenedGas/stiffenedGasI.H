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

#include "stiffenedGas.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#if (OPENFOAM >= 1912)
#define small SMALL
#endif

template<class Specie>
inline Foam::stiffenedGas<Specie>::stiffenedGas
(
    const Specie& sp,
    const scalar& pInf
)
:
    Specie(sp),
    pInf_(pInf)
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::stiffenedGas<Specie>::stiffenedGas
(
    const word& name,
    const stiffenedGas& pg
)
:
    Specie(name, pg),
    pInf_(pg.pInf_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::stiffenedGas <Specie>>
Foam::stiffenedGas<Specie>::clone() const
{
    return autoPtr<stiffenedGas<Specie>>
    (
        new stiffenedGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::stiffenedGas<Specie>>
Foam::stiffenedGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<stiffenedGas<Specie>>
    (
        new stiffenedGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::rho(scalar p, scalar T) const
{
    return (p + pInf_)/(this->R()*T);
}



template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::H(scalar p, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::Cp(scalar p, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::E(const scalar p, const scalar T) const
{
    return pInf_/this->rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::Cv(scalar p, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::S(scalar p, scalar T) const
{
    return -this->R()*log((p + pInf_)/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::psi(scalar p, scalar T) const
{
    const scalar Z = this->Z(p, T);
    return 1.0/(Z*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::Z(scalar p, scalar T) const
{
    return p/(p + pInf_);
}


template<class Specie>
inline Foam::scalar Foam::stiffenedGas<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::stiffenedGas<Specie>::operator+=
(
    const stiffenedGas<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();

        pInf_ = Y1*pInf_ + Y2*pg.pInf_;
    }
}


template<class Specie>
inline void Foam::stiffenedGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::stiffenedGas<Specie> Foam::operator+
(
    const stiffenedGas<Specie>& pg1,
    const stiffenedGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return stiffenedGas<Specie>
        (
            sp,
            pg1.pInf_
        );
    }
    else
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();
        
        return stiffenedGas<Specie>
        (
            sp,
            Y1*pg1.pInf_ + Y2*pg2.pInf_
        );
    }
}


template<class Specie>
Foam::stiffenedGas<Specie> Foam::operator*
(
    const scalar s,
    const stiffenedGas<Specie>& pg
)
{
    return stiffenedGas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.pInf_
    );
}


template<class Specie>
Foam::stiffenedGas<Specie> Foam::operator==
(
    const stiffenedGas<Specie>& pg1,
    const stiffenedGas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();

    const scalar pInf = Y2*pg2.pInf_ - Y1*pg1.pInf_;

    return stiffenedGas<Specie>
    (
        sp,
        pInf
    );
}


// ************************************************************************* //
