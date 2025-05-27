/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::sutherlandPolynomialTransport<Thermo, PolySize>::calcCoeffs
(
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2
)
{
    scalar rootT1 = sqrt(T1);
    scalar mu1rootT2 = mu1*sqrt(T2);
    scalar mu2rootT1 = mu2*rootT1;

    Ts_ = (mu2rootT1 - mu1rootT2)/(mu1rootT2/T1 - mu2rootT1/T2);

    As_ = mu1*(1.0 + Ts_/T1)/rootT1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::sutherlandPolynomialTransport<Thermo, PolySize>::sutherlandPolynomialTransport
(
    const Thermo& t,
    const scalar As,
    const scalar Ts,
    const Polynomial<PolySize>& kappaCoeffs
)
:
    Thermo(t),
    As_(As),
    Ts_(Ts),
    kappaCoeffs_(kappaCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::sutherlandPolynomialTransport<Thermo, PolySize>::sutherlandPolynomialTransport
(
    const Thermo& t,
    const scalar mu1, const scalar T1,
    const scalar mu2, const scalar T2,
    const Polynomial<PolySize>& kappaCoeffs
)
:
    Thermo(t),
    kappaCoeffs_(kappaCoeffs)
{
    calcCoeffs(mu1, T1, mu2, T2);
}


template<class Thermo, int PolySize>
inline Foam::sutherlandPolynomialTransport<Thermo, PolySize>::sutherlandPolynomialTransport
(
    const word& name,
    const sutherlandPolynomialTransport& stp
)
:
    Thermo(name, stp),
    As_(stp.As_),
    Ts_(stp.Ts_),
    kappaCoeffs_(stp.kappaCoeffs)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::sutherlandPolynomialTransport<Thermo, PolySize>>
Foam::sutherlandPolynomialTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<sutherlandPolynomialTransport<Thermo, PolySize>>::New(*this);
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::sutherlandPolynomialTransport<Thermo, PolySize>>
Foam::sutherlandPolynomialTransport<Thermo, PolySize>::New
(
    const dictionary& dict
)
{
    return autoPtr<sutherlandPolynomialTransport<Thermo, PolySize>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::sutherlandPolynomialTransport<Thermo, PolySize>::mu
(
    const scalar p,
    const scalar T
) const
{
    return As_*::sqrt(T)/(1.0 + Ts_/T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::sutherlandPolynomialTransport<Thermo, PolySize>::kappa
(
    const scalar p, const scalar T
) const
{
    return kappaCoeffs_.value(T);
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::sutherlandPolynomialTransport<Thermo, PolySize>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline void Foam::sutherlandPolynomialTransport<Thermo, PolySize>::operator+=
(
    const sutherlandPolynomialTransport<Thermo, PolySize>& stp
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(stp);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = stp.Y()/this->Y();

        As_ = Y1*As_ + Y2*stp.As_;
        Ts_ = Y1*Ts_ + Y2*stp.Ts_;

        kappaCoeffs_ = Y1*kappaCoeffs_ + Y2*stp.kappaCoeffs_;
    }
}


template<class Thermo, int PolySize>
inline void Foam::sutherlandPolynomialTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::sutherlandPolynomialTransport<Thermo, PolySize> Foam::operator+
(
    const sutherlandPolynomialTransport<Thermo, PolySize>& stp1,
    const sutherlandPolynomialTransport<Thermo, PolySize>& stp2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(stp1) + static_cast<const Thermo&>(stp2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return sutherlandPolynomialTransport<Thermo, PolySize>
        (
            t,
            stp1.As_,
            stp1.Ts_,
            stp1.kappaCoeffs_
        );
    }
    else
    {
        scalar Y1 = stp1.Y()/t.Y();
        scalar Y2 = stp2.Y()/t.Y();

        return sutherlandPolynomialTransport<Thermo, PolySize>
        (
            t,
            Y1*stp1.As_ + Y2*stp2.As_,
            Y1*stp1.Ts_ + Y2*stp2.Ts_,
            Y1*stp1.kappaCoeffs_ + Y2*stp2.kappaCoeffs_
        );
    }
}


template<class Thermo, int PolySize>
inline Foam::sutherlandPolynomialTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const sutherlandPolynomialTransport<Thermo, PolySize>& stp
)
{
    return sutherlandPolynomialTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(stp),
        stp.As_,
        stp.Ts_,
        stp.kappaCoeffs_
    );
}


// ************************************************************************* //
