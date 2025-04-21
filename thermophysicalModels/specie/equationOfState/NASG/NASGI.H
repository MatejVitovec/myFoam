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

#include "NASG.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#if (OPENFOAM >= 1912)
#define small SMALL
#endif

template<class Specie>
inline Foam::NASG<Specie>::NASG
(
    const Specie& sp,
    const scalar& pInf,
    const scalar& b
)
:
    Specie(sp),
    pInf_(pInf),
    b_(b)
{
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::NASG<Specie>::NASG
(
    const word& name,
    const NASG& pg
)
:
    Specie(name, pg),
    pInf_(pg.pInf_),
    b_(pg.b_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::NASG <Specie>>
Foam::NASG<Specie>::clone() const
{
    return autoPtr<NASG<Specie>>
    (
        new NASG<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::NASG<Specie>>
Foam::NASG<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<NASG<Specie>>
    (
        new NASG<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::rho(scalar p, scalar T) const
{
    return (p + pInf_)/(this->R()*T + (p + pInf_)*b_);
}



template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::H(scalar p, scalar T) const
{
    return p*b_;
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::Cp(scalar p, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::E(const scalar p, const scalar T) const
{
    return (pInf_*this->R()*T)/(p + pInf_);
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::Cv(scalar p, scalar T) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::S(scalar p, scalar T) const
{
    //return -this->R()*log((p + pInf_)/(Pstd + pInf_));
    return -this->R()*log((p + pInf_)/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::psi(scalar p, scalar T) const
{
    //const scalar Z = this->Z(p, T);
    //return 1.0/(Z*this->R()*T);
    return rho(p, T)/p;
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::Z(scalar p, scalar T) const
{
    return (p*b_)/(this->R()*T) - pInf_/(p + pInf_);
}


template<class Specie>
inline Foam::scalar Foam::NASG<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::NASG<Specie>::operator+=
(
    const NASG<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();

        pInf_ = Y1*pInf_ + Y2*pg.pInf_;
        b_ = Y1*b_ + Y2*pg.b_;
    }
}


template<class Specie>
inline void Foam::NASG<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::NASG<Specie> Foam::operator+
(
    const NASG<Specie>& pg1,
    const NASG<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < small)
    {
        return NASG<Specie>
        (
            sp,
            pg1.pInf_,
            pg1.b_
        );
    }
    else
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();
        
        return NASG<Specie>
        (
            sp,
            Y1*pg1.pInf_ + Y2*pg2.pInf_,
            Y1*pg1.b_ + Y2*pg2.b_
        );
    }
}


template<class Specie>
Foam::NASG<Specie> Foam::operator*
(
    const scalar s,
    const NASG<Specie>& pg
)
{
    return NASG<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.pInf_,
        pg.b_
    );
}


template<class Specie>
Foam::NASG<Specie> Foam::operator==
(
    const NASG<Specie>& pg1,
    const NASG<Specie>& pg2
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
    const scalar b = Y2*pg2.b_ - Y1*pg1.b_;

    return NASG<Specie>
    (
        sp,
        pInf,
        b
    );
}


// ************************************************************************* //
