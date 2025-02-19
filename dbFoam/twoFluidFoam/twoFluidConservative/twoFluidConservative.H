/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    dbFoam

Description
    Density-based compressible steady-state & transient flow solver.

Author
    Matej Vitovec

\*---------------------------------------------------------------------------*/

#ifndef twoFluidConservative_H
#define twoFluidConservative_H

#include "fvCFD.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class twoFluidConservative Declaration
\*---------------------------------------------------------------------------*/

class twoFluidConservative
{        
    volScalarField alphaRho1_;
    volScalarField alphaRho2_;
    volVectorField alphaRhoU1_;
    volVectorField alphaRhoU2_;
    volScalarField epsilon1_;
    volScalarField epsilon2_;

    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    twoFluidConservative(const twoFluidConservative&);
        
    //- Disallow default bitwise assignment
    void operator=(const twoFluidConservative&);
    

public:

    // Constructors

    //- Construct from components
    twoFluidConservative
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
    );


    //- Destructor
    virtual ~twoFluidConservative()
    {}


    // Memeber functions


    const volScalarField& alphaRho1() const
    {
        return alphaRho1_;
    }

    volScalarField& alphaRho1()
    {
        return alphaRho1_;
    }

    const volScalarField& alphaRho2() const
    {
        return alphaRho2_;
    }

    volScalarField& alphaRho2()
    {
        return alphaRho2_;
    }

    const volVectorField& alphaRhoU1() const
    {
        return alphaRhoU1_;
    }

    volVectorField& alphaRhoU1()
    {
        return alphaRhoU1_;
    }

    const volVectorField& alphaRhoU2() const
    {
        return alphaRhoU2_;
    }

    volVectorField& alphaRhoU2()
    {
        return alphaRhoU2_;
    }

    const volScalarField& epsilon1() const
    {
        return epsilon1_;
    }

    volScalarField& epsilon1()
    {
        return epsilon1_;
    }

    const volScalarField& epsilon2() const
    {
        return epsilon2_;
    }

    volScalarField& epsilon2()
    {
        return epsilon2_;
    }

    void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //