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

#include "sutherlandPolynomialTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::scalar Foam::sutherlandPolynomialTransport<Thermo, PolySize>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").get<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::sutherlandPolynomialTransport<Thermo, PolySize>::sutherlandPolynomialTransport(const dictionary& dict)
:
    Thermo(dict),
    As_(readCoeff("As", dict)),
    Ts_(readCoeff("Ts", dict)),
    kappaCoeffs_(dict.subDict("transport").lookup(coeffsName("kappa")))
{}


template<class Thermo, int PolySize>
Foam::sutherlandPolynomialTransport<Thermo, PolySize>::sutherlandPolynomialTransport
(
    const Thermo& t,
    const dictionary& dict
)
:
    Thermo(t),
    As_(readCoeff("As", dict)),
    Ts_(readCoeff("Ts", dict)),
    kappaCoeffs_(dict.subDict("transport").lookup(coeffsName("kappa")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
void Foam::sutherlandPolynomialTransport<Thermo, PolySize>::write(Ostream& os) const
{
    os.beginBlock(this->specie::name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.writeEntry("As", As_);
        os.writeEntry("Ts", Ts_);
        os.writeEntry(coeffsName("kappa"), kappaCoeffs_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const sutherlandPolynomialTransport<Thermo, PolySize>& stp
)
{
    stp.write(os);
    return os;
}


// ************************************************************************* //
