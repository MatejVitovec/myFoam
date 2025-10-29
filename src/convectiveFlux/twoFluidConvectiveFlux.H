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

#ifndef twoFluidConvectiveFlux_H
#define twoFluidConvectiveFlux_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "riemannSolver.H"
#include "rhoThermo.H"
#include "gasProperties.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class twoFluidConvectiveFlux Declaration
\*---------------------------------------------------------------------------*/

class twoFluidConvectiveFlux
{

    autoPtr<riemannSolver> fluxSolver_;
        
    const fvMesh& mesh_;
        
    const volScalarField& p_;
    const volScalarField& alpha_;
    const volVectorField& U1_;
    const volVectorField& U2_;
    const volScalarField& T1_;
    const volScalarField& T2_;

    //gasProperties& gasProps1_;
    //gasProperties& gasProps1_;
    autoPtr<gasProperties> gasProps1_;
    autoPtr<gasProperties> gasProps2_;

    surfaceScalarField alpha_pos_;
    surfaceScalarField alpha_neg_;

    surfaceScalarField alphaRhoFlux1_pos_;
    surfaceScalarField alphaRhoFlux1_neg_;
    surfaceScalarField alphaRhoFlux2_pos_;
    surfaceScalarField alphaRhoFlux2_neg_;

    surfaceVectorField alphaRhoUFlux1_pos_;
    surfaceVectorField alphaRhoUFlux1_neg_;
    surfaceVectorField alphaRhoUFlux2_pos_;
    surfaceVectorField alphaRhoUFlux2_neg_;

    surfaceScalarField alphaRhoEFlux1_pos_;
    surfaceScalarField alphaRhoEFlux1_neg_;
    surfaceScalarField alphaRhoEFlux2_pos_;
    surfaceScalarField alphaRhoEFlux2_neg_;
    
    
    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    twoFluidConvectiveFlux(const twoFluidConvectiveFlux&);
        
    //- Disallow default bitwise assignment
    void operator=(const twoFluidConvectiveFlux&);

    const gasProperties& gasProps1() const
    {
        return *gasProps1_;
    }

    const gasProperties& gasProps2() const
    {
        return *gasProps2_;
    }
    

public:

    // Constructors

    //- Construct from components
    twoFluidConvectiveFlux
    (
        const volScalarField& p,
        const volScalarField& alpha,
        const volVectorField& U1,
        const volVectorField& U2,
        const volScalarField& T1,
        const volScalarField& T2,
        const rhoThermo& thermo1,
        const rhoThermo& thermo2
    );


    //- Destructor
    virtual ~twoFluidConvectiveFlux()
    {}


    // Memeber functions

    const fvMesh& mesh() const
    {
        return mesh_;
    }

    const surfaceScalarField& alpha_pos() const
    {
        return alpha_pos_;
    }

    const surfaceScalarField& alpha_neg() const
    {
        return alpha_neg_;
    }

    const tmp<surfaceScalarField> alpha1_pos() const
    {
        return 1.0 - alpha_pos_;
    }

    const tmp<surfaceScalarField> alpha1_neg() const
    {
        return 1.0 - alpha_neg_;
    }

    const surfaceScalarField& alpha2_pos() const
    {
        return alpha_pos_;
    }

    const surfaceScalarField& alpha2_neg() const
    {
        return alpha_neg_;
    }

    const surfaceScalarField& alphaRhoFlux1_pos() const
    {
        return alphaRhoFlux1_pos_;
    }

    const surfaceScalarField& alphaRhoFlux1_neg() const
    {
        return alphaRhoFlux1_neg_;
    }

    const surfaceScalarField& alphaRhoFlux2_pos() const
    {
        return alphaRhoFlux2_pos_;
    }

    const surfaceScalarField& alphaRhoFlux2_neg() const
    {
        return alphaRhoFlux2_neg_;
    }

    const surfaceVectorField& alphaRhoUFlux1_pos() const
    {
        return alphaRhoUFlux1_pos_;
    }

    const surfaceVectorField& alphaRhoUFlux1_neg() const
    {
        return alphaRhoUFlux1_neg_;
    }

    const surfaceVectorField& alphaRhoUFlux2_pos() const
    {
        return alphaRhoUFlux2_pos_;
    }

    const surfaceVectorField& alphaRhoUFlux2_neg() const
    {
        return alphaRhoUFlux2_neg_;
    }

    const surfaceScalarField& alphaRhoEFlux1_pos() const
    {
        return alphaRhoEFlux1_pos_;
    }

    const surfaceScalarField& alphaRhoEFlux1_neg() const
    {
        return alphaRhoEFlux1_neg_;
    }

    const surfaceScalarField& alphaRhoEFlux2_pos() const
    {
        return alphaRhoEFlux2_pos_;
    }

    const surfaceScalarField& alphaRhoEFlux2_neg() const
    {
        return alphaRhoEFlux2_neg_;
    }

    void computeFlux();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //