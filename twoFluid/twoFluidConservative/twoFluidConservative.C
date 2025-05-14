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

#include "fvCFD.H"
#include "twoFluidConservative.H"
#include <iostream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::twoFluidConservative::twoFluidConservative
(
    const volScalarField& p,
    const volScalarField& alpha,
    const volVectorField& U1,
    const volVectorField& U2,
    const volScalarField& T1,
    const volScalarField& T2,
    const rhoThermo& thermo1,
    const rhoThermo& thermo2,
    const volScalarField& pInt
)
:
    alphaRho1_
    (
        IOobject
        (
            "alphaRho1",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (1.0 - alpha)*thermo1.rho()
    ),
    alphaRho2_
    (
        IOobject
        (
            "alphaRho2",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha*thermo2.rho()
    ),
    alphaRhoU1_
    (
        IOobject
        (
            "alphaRhoU1",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRho1_*U1
    ),
    alphaRhoU2_
    (
        IOobject
        (
            "alphaRhoU2",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alphaRho2_*U2
    ),
    epsilon1_
    (
        IOobject
        (
            "epsilon1",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (1.0 - alpha)*(thermo1.rho()*(thermo1.he() + 0.5*Foam::magSqr(U1)) + pInt)
    ),
    epsilon2_
    (
        IOobject
        (
            "epsilon2",
            p.time().timeName(),
            p.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha*(thermo2.rho()*(thermo2.he() + 0.5*Foam::magSqr(U2)) + pInt)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //