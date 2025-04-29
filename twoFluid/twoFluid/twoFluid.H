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

#include "fvCFD.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "rhoThermo.H"
#include "gasProperties.H"

#include "twoFluidConservative.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace TwoFluidFoam
{

/*---------------------------------------------------------------------------*\
                        Class twoFluid Declaration
\*---------------------------------------------------------------------------*/

class twoFluid
:
    public IOdictionary
{
    const fvMesh& mesh_;

    const scalar epsilon_;
    const scalar epsilonMin_;
    const scalar epsilonMax_;

    word fluid1Name_;
    word fluid2Name_;

    autoPtr<rhoThermo> pthermo1_;
    autoPtr<rhoThermo> pthermo2_;
    rhoThermo& thermo1_;
    rhoThermo& thermo2_;

    autoPtr<gasProperties> pGasProps1_;
    autoPtr<gasProperties> pGasProps2_;
    gasProperties& gasProps1_;
    gasProperties& gasProps2_;

    volScalarField& p_;
    volScalarField  alpha_;
    volScalarField  alpha2_;
    volVectorField  U1_;
    volVectorField  U2_;
    volScalarField& T1_;
    volScalarField& T2_;

    volScalarField a1_;
    volScalarField a2_;

    volScalarField pInt_;
    
    twoFluidConservative conservative_;
            
    
    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    twoFluid(const twoFluid&);
        
    //- Disallow default bitwise assignment
    void operator=(const twoFluid&);

    const gasProperties& gas1() const
    {
        return gasProps1_;
    }

    const gasProperties& gas2() const
    {
        return gasProps2_;
    }

    /*void primitiveFromConservative
    (
        scalar& p,
        scalar& alpha,
        vector& U1,
        vector& U2,
        scalar& T1,
        scalar& T2,
        const scalar alphaRho1,
        const scalar alphaRho2,
        const vector alphaRhoU1,
        const vector alphaRhoU2,
        const scalar epsilon1,
        const scalar epsilon2,
        const scalar pIntOld
    );*/

    void correctSoundSpeeds();
    
public:

    // Constructors

    //- Construct from components
    /*twoFluid
    (
        const volScalarField& p,
        const volScalarField& alpha,
        const volVectorField& U1,
        const volVectorField& U2,
        const volScalarField& T1,
        const volScalarField& T2,
        rhoThermo& thermo1,
        rhoThermo& thermo2
    );*/

    twoFluid(const fvMesh& mesh);


    //- Destructor
    virtual ~twoFluid()
    {}

    void primitiveFromConservative
    (
        scalar& p,
        scalar& alpha,
        vector& U1,
        vector& U2,
        scalar& T1,
        scalar& T2,
        const scalar alphaRho1,
        const scalar alphaRho2,
        const vector alphaRhoU1,
        const vector alphaRhoU2,
        const scalar epsilon1,
        const scalar epsilon2,
        const scalar pIntOld
    ) const;

    void blendVanishingFluid
    (
        scalar& alpha,
        vector& U1,
        vector& U2,
        scalar& T1,
        scalar& T2
    ) const;


    // Memeber functions

    scalar epsilon() const
    {
        return epsilon_;
    } 

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

    const volScalarField& alpha2() const
    {
        return alpha2_;
    }

    volScalarField& alpha2() 
    {
        return alpha2_; //TODO predelat v condenzaci a pak odstranit
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

    const volScalarField& a1() const
    {
        return a1_;
    }

    const volScalarField& a2() const
    {
        return a2_;
    }

    const volScalarField& pInt() const
    {
        return pInt_;
    }

    volScalarField& pInt()
    {
        return pInt_;
    }

    const rhoThermo& thermo1() const
    {
        return thermo1_;
    }

    rhoThermo& thermo1()
    {
        return thermo1_;
    }

    const rhoThermo& thermo2() const
    {
        return thermo2_;
    }

    rhoThermo& thermo2()
    {
        return thermo2_;
    }

    const gasProperties& gasProps1() const
    {
        return gasProps1_;
    }

    const gasProperties& gasProps2() const
    {
        return gasProps2_;
    }

    const twoFluidConservative& conservative() const
    {
        return conservative_;
    }

    twoFluidConservative& conservative()
    {
        return conservative_;
    }

    void correct();

    void blendVanishingFluid();

    void correctBoundaryCondition();

    void correctThermo();

    void correctInterfacialPressure();

    void correctConservative();

    tmp<surfaceScalarField> amaxSf() const;
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace TwoFluidFoam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //