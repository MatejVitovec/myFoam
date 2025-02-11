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

#ifndef twoFluid_H
#define twoFluid_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "riemannSolver.H"
#include "rhoThermo.H"
#include "gasProperties.H"

#include "twoFluidConservative.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace twoFluidFoam
{

/*---------------------------------------------------------------------------*\
                        Class twoFluid Declaration
\*---------------------------------------------------------------------------*/

class twoFluid
{        
    const fvMesh& mesh_;
        
    volScalarField& p_;
    volScalarField& alpha_;
    volVectorField& U1_;
    volVectorField& U2_;
    volScalarField& T1_;
    volScalarField& T2_;
        
    rhoThermo& thermo1_;
    rhoThermo& thermo2_;
    autoPtr<gasProperties> gasProps1_;
    autoPtr<gasProperties> gasProps2_;

    volScalarField& pInt_;

    twoFluidConservative& conservative_;
            
    
    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    twoFluid(const twoFluid&);
        
    //- Disallow default bitwise assignment
    void operator=(const twoFluid&);

    const gasProperties& gas1() const
    {
        return *gasProps1_;
    }

    const gasProperties& gas2() const
    {
        return *gasProps2_;
    }

    
public:

    // Constructors

    //- Construct from components
    twoFluid
    (
        const volScalarField& p,
        const volScalarField& alpha,
        const volVectorField& U1,
        const volVectorField& U2,
        const volScalarField& T1,
        const volScalarField& T2,
        rhoThermo& thermo1,
        rhoThermo& thermo2
    );

    twoFluid(const fvMesh& mesh);


    //- Destructor
    virtual ~twoFluid()
    {}


    // Memeber functions

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    const volScalarField& p() const
    {
        return p_;
    }

    volScalarField& p()
    {
        return p_;
    }

    const volScalarField& alpha() const
    {
        return alpha_;
    }

    volScalarField& alpha()
    {
        return alpha_;
    }

    const volVectorField& U1() const
    {
        return U1_;
    }

    volVectorField& U1()
    {
        return U1_;
    }

    const volVectorField& U2() const
    {
        return U2_;
    }

    volVectorField& U2()
    {
        return U2_;
    }

    const volScalarField& T1() const
    {
        return T1_;
    }

    volScalarField& T1()
    {
        return T1_;
    }

    const volScalarField& T2() const
    {
        return T2_;
    }

    volScalarField& T2()
    {
        return T2_;
    }

    void correct();

    void blend();

    void correctBoundaryCondition();

    void correctThermo();

    void correctInterfacialPressure();

    void correctConservative();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoFluidFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //