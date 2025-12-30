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

#ifndef compressibleMixture_H
#define compressibleMixture_H

#include "fvCFD.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "rhoThermo.H"
#include "gasProperties.H"

#include "compressibleMixtureConservative.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class compressibleMixture Declaration
\*---------------------------------------------------------------------------*/

class compressibleMixture
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
    volVectorField  U_;
    volScalarField& T_;
    volScalarField  w_;

    volScalarField a_;
    
    compressibleMixtureConservative conservative_;
            
    
    // Private Member Functions
    
    //- Disallow default bitwise copy construct
    compressibleMixture(const compressibleMixture&);
        
    //- Disallow default bitwise assignment
    void operator=(const compressibleMixture&);

    const gasProperties& gas1() const
    {
        return gasProps1_;
    }

    const gasProperties& gas2() const
    {
        return gasProps2_;
    }


    void correctSoundSpeeds();
    
public:

    // Constructors

    //- Construct from components
    compressibleMixture
    (
        volScalarField& p,
        volVectorField& U,
        volScalarField& T,
        volScalarField& w,
        rhoThermo& thermo1,
        rhoThermo& thermo2
    );

    compressibleMixture(const fvMesh& mesh);


    //- Destructor
    virtual ~compressibleMixture()
    {}

    void primitiveFromConservative
    (
        scalar& p,
        vector& U,
        scalar& T,
        scalar& w,
        const scalar rho,
        const vector rhoU,
        const scalar rhoE,
        const scalar rhoW
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

    const volScalarField& w() const
    {
        return w_;
    }

    volScalarField& w()
    {
        return w_;
    }

    const volVectorField& U() const
    {
        return U_;
    }

    volVectorField& U()
    {
        return U_;
    }

    const volScalarField& T() const
    {
        return T_;
    }

    volScalarField& T()
    {
        return T_;
    }

    const volScalarField& a() const
    {
        return a_;
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

    tmp<volScalarField> rho() const
    {
        //return 1.0/(w_/thermo1_.rho() + (1.0 - w)/thermo2_.rho());
        return thermo1_.rho()/(1.0 - w_);
    }

    tmp<volScalarField> alpha() const
    {
        //return thermo1_.rho()/(thermo1_rho()*w_ + (1.0 - w_)*thermo2_rho());
        return w_*rho()/thermo2_.rho();
    }

    volScalarField& rho1() const
    {
        return thermo1_.rho();
    }

    volScalarField& rho2() const
    {
        return thermo2_.rho();
    }

    volScalarField& he1() const
    {
        return thermo1_.he();
    }

    volScalarField& he2() const
    {
        return thermo2_.he();
    }

    const gasProperties& gasProps1() const
    {
        return gasProps1_;
    }

    const gasProperties& gasProps2() const
    {
        return gasProps2_;
    }

    const compressibleMixtureConservative& conservative() const
    {
        return conservative_;
    }

    compressibleMixtureConservative& conservative()
    {
        return conservative_;
    }

    void correct();

    void correctBoundaryCondition();

    void correctThermo();

    void correctConservative();

    tmp<surfaceScalarField> amaxSf() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //