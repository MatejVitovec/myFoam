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

#ifndef compressibleMixtureConservative_H
#define compressibleMixtureConservative_H

#include "fvCFD.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class compressibleMixtureConservative Declaration
\*---------------------------------------------------------------------------*/

class compressibleMixtureConservative
{        
    volScalarField rho_;
    volVectorField rhoU_;
    volScalarField rhoE_;
    volScalarField rhow_;

    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    compressibleMixtureConservative(const compressibleMixtureConservative&);
        
    //- Disallow default bitwise assignment
    void operator=(const compressibleMixtureConservative&);
    

public:

    // Constructors

    //- Construct from components
    compressibleMixtureConservative
    (
        const volScalarField& p,
        const volVectorField& U,
        const volScalarField& T,
        const volScalarField& w,
        const rhoThermo& thermo1,
        const rhoThermo& thermo2
    );


    //- Destructor
    virtual ~compressibleMixtureConservative()
    {}


    // Memeber functions

    const volScalarField& rho() const
    {
        return rho_;
    }

    volScalarField& rho()
    {
        return rho_;
    }

    const volVectorField& rhoU() const
    {
        return rhoU_;
    }

    volVectorField& rhoU()
    {
        return rhoU_;
    }

    const volScalarField& rhoE() const
    {
        return rhoE_;
    }

    volScalarField& rhoE()
    {
        return rhoE_;
    }

    const volScalarField& rhow() const
    {
        return rhow_;
    }

    volScalarField& rhow()
    {
        return rhow_;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //