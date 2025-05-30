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

#ifndef convectiveFlux_H
#define convectiveFlux_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "riemannSolver.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class convectiveFlux Declaration
\*---------------------------------------------------------------------------*/

class convectiveFlux
{
    // Private data

        autoPtr<riemannSolver> fluxSolver_;
        
        //- Reference to mesh
        const fvMesh& mesh_;
        
        // Reference to primitive fields
        
            //- Static pressure
            const volScalarField& p_;
            
            //- Velocity
            const volVectorField& U_;
            
            //- Static temperature
            const volScalarField& T_;
        
        //- Reference to the thermophysicalModel
        psiThermo& thermo_;
    

        // Fluxes
        
            //- Density flux
            surfaceScalarField rhoFlux_;
            
            //- Velocity flux
            surfaceVectorField rhoUFlux_;
            
            //- Energy flux
            surfaceScalarField rhoEFlux_;
    
    // Private Member Functions
    
        //- Disallow default bitwise copy construct
        convectiveFlux(const convectiveFlux&);
        
        //- Disallow default bitwise assignment
        void operator=(const convectiveFlux&);
        
        
        //- Return internal field of mass flux
        const scalarField& rhoFluxI() const
        {
            return rhoFlux_.internalField();
        }
    
    
        //- Return internal field of momentum flux
        const vectorField& rhoUFluxI() const
        {
            return rhoUFlux_.internalField();
        }
    

        //- Return access to internal field of energy flux
        const scalarField& rhoEFluxI() const
        {
            return rhoEFlux_.internalField();
        }
    

public:

    // Constructors

        //- Construct from components
        convectiveFlux
        (
            const volScalarField& p,
            const volVectorField& U,
            const volScalarField& T,
            psiThermo& thermo
        );


    //- Destructor
    virtual ~convectiveFlux()
    {}


    // Memeber functions

        const fvMesh& mesh() const
        {
            return mesh_;
        }

        // Return fluxes

            //- Return continuity equation flux
            const surfaceScalarField& rhoFlux() const
            {
                return rhoFlux_;
            }

            //- Return momentum equation flux
            const surfaceVectorField& rhoUFlux() const
            {
                return rhoUFlux_;
            }

            //- Return energy equation flux
            const surfaceScalarField& rhoEFlux() const
            {
                return rhoEFlux_;
            }


        // Return residuals

            //- Return mass residual
            tmp<scalarField> rhoResidual() const
            {
                return fvc::div(rhoFlux_)().internalField();
            }

            //- Return momentum residual
            tmp<vectorField> rhoUResidual() const
            {
                return fvc::div(rhoUFlux_)().internalField();
            }

            //- Return energy residual
            tmp<scalarField> rhoEResidual() const
            {
                return fvc::div(rhoEFlux_)().internalField();
            }


        //- Update fluxes based on current state
        void computeFlux();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //