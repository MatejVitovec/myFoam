/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvCFD.H"
#include "fvMesh.H"
#include "rhoThermo.H"

#include "twoFluidFvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TwoFluidFoam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void primitiveFromConservative
(
    volScalarField& p,
    volScalarField& alpha,
    volVectorField& U1,
    volVectorField& U2,
    volScalarField& T1,
    volScalarField& T2,
    const volScalarField& alphaRho1,
    const volScalarField& alphaRho2,
    const volVectorField& alphaRhoU1,
    const volVectorField& alphaRhoU2,
    const volScalarField& alphaRhoE1,
    const volScalarField& alphaRhoE2,
    const rhoThermo& thermo1,
    const rhoThermo& thermo2,
    const fvMesh& mesh
)
{
    const auto& cells = mesh.C();

    forAll(cells, celli)
    {
        //TODO
    }
}



void updateIntPressure
(
    volScalarField& pInt,
    const volScalarField& alpha,
    const volVectorField& U1,
    const volVectorField& U2,
    const rhoThermo& thermo1,
    const rhoThermo& thermo2
)
{
    const scalar sigma = 2.0;
    const scalar epsilonP = 0.01;

    const auto& rho1 = thermo1.rho();
    const auto& rho2 = thermo2.rho();

    pInt = sigma*((alpha*rho1*(1.0 - alpha)*rho2)/
        (alpha*rho2 + (1.0 - alpha)*rho1))*magSqr(U2 - U1);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TwoFluidFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
