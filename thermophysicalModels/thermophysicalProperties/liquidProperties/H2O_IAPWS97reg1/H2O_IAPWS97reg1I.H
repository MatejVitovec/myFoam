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

\*---------------------------------------------------------------------------*/

Foam::scalar Foam::H2O_IAPWS97reg1::IAPWS_sigma::f(scalar, scalar T) const
{
    const scalar Tc  = 647.096;
    const scalar tau = max(1.0 - T/Tc, 0.0);

    return 235.8e-3*pow(tau, 1.256)*(1 - 0.625*tau);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::rho(scalar p, scalar T) const
{
    return 1.0/IF97::reg1::v(p, T);
}

inline Foam::scalar Foam::H2O_IAPWS97reg1::psi(scalar p, scalar T) const
{
    return p/(rho(p, T)*T*(Foam::constant::thermodynamic::RR/W()));
}

inline Foam::scalar Foam::H2O_IAPWS97reg1::CpMCv(scalar p, scalar T) const
{
    return IF97::reg1::cp(p, T) - IF97::reg1::cv(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::pv(scalar p, scalar T) const
{
    return pv_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::hl(scalar p, scalar T) const
{
    return hl_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::Cp(scalar p, scalar T) const
{
    return IF97::reg1::cp(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::h(scalar p, scalar T) const
{
    return IF97::reg1::h(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::Cpg(scalar p, scalar T) const
{
    return IF97::reg2meta::cp(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::B(scalar p, scalar T) const
{
    return B_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::mu(scalar p, scalar T) const
{
    return mu_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::mug(scalar p, scalar T) const
{
    return IF97::mu(1/IF97::reg2meta::v(p, T), T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::kappa(scalar p, scalar T) const
{
    return kappa_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::kappag(scalar p, scalar T) const
{
    return IF97::lambda(1/IF97::reg2meta::v(p, T), T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::sigma(scalar p, scalar T) const
{
    return sigma_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::D(scalar p, scalar T) const
{
    return D_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS97reg1::D(scalar p, scalar T, scalar Wb) const
{
    return D_.f(p, T, Wb);
}


// ************************************************************************* //
