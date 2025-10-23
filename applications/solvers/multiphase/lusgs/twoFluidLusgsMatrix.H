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

#ifndef twoFluidLusgsMatrix_H
#define twoFluidLusgsMatrix_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "riemannSolver.H"
#include "psiThermo.H"
#include "gasProperties.H"

#include <eigen3/Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class twoFluidLusgsMatrix Declaration
\*---------------------------------------------------------------------------*/

class twoFluidLusgsMatrix
{
    // Private data

        PartialPivLU<Matrix<double,10,10>> D;
        
        //- Reference to mesh
        const fvMesh& mesh_;
        
        // Reference to primitive fields
        
            //- Static pressure
            const volScalarField& p_;
            
            //- Velocity
            const volVectorField& U_;
            
            //- Static temperature
            const volScalarField& T_;
    


    
    // Private Member Functions
    
        //- Disallow default bitwise copy construct
        twoFluidLusgsMatrix(const twoFluidLusgsMatrix&);
        
        //- Disallow default bitwise assignment
        void operator=(const twoFluidLusgsMatrix&);

public:

    // Constructors

        //- Construct from components
        twoFluidLusgsMatrix
        (
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& T,
            psiThermo& thermo,
            gasProperties& gasProps
        );


    //- Destructor
    virtual ~twoFluidLusgsMatrix()
    {}


    // Memeber functions

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Update fluxes based on current state
        void computeFlux();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //