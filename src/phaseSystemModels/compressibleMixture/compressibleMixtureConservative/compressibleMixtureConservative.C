/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "compressibleMixtureConservative.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::compressibleMixtureConservative::compressibleMixtureConservative
(
    const volScalarField& p,
    const volVectorField& U,
    const volScalarField& T,
    const volScalarField& w,
    const rhoThermo& thermo1,
    const rhoThermo& thermo2
)
:
    rho_
    (
        IOobject
        (
            "rho",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        1.0/((1.0 - w)/thermo1.rho() + w/thermo2.rho()) //TODO
    ),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*U
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*((1.0 - w)*thermo1.he() + w*thermo2.he() + 0.5*magSqr(U))
    ),
    rhow_
    (
        IOobject
        (
            "rhow",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho_*w
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //