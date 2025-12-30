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
#include "compressibleMixture.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::compressibleMixture::compressibleMixture
(
    volScalarField& p,
    volVectorField& U,
    volScalarField& T,
    volScalarField& w,
    rhoThermo& thermo1,
    rhoThermo& thermo2
)
:
    IOdictionary
    (
        IOobject
        (
            "compressibleMixtureProperties",
            p.mesh().time().constant(),
            p.mesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(p.mesh()),
    epsilon_(this->lookupOrDefault<scalar>("epsilon", 1.0e-10)),
    epsilonMin_(this->lookupOrDefault<scalar>("epsilonMin", epsilon_*0.1)),
    epsilonMax_(this->lookupOrDefault<scalar>("epsilonMax", epsilon_*100)),
    pthermo1_(&thermo1),
    pthermo2_(&thermo2),
    thermo1_(pthermo1_()),
    thermo2_(pthermo2_()),
    pGasProps1_(gasProperties::New(thermo1_)),
    pGasProps2_(gasProperties::New(thermo2_)),
    gasProps1_(pGasProps1_()),
    gasProps2_(pGasProps2_()),
    p_(p),
    U_(U),
    T_(T),
    w_(w),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimLength/dimTime
    ),
    conservative_
    (
        p_,
        U_,
        T_,
        w_,
        thermo1_,
        thermo2_
    )
{
    correctBoundaryCondition();
    correctThermo();
    correctConservative();
}


Foam::compressibleMixture::compressibleMixture(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "compressibleMixtureProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    epsilon_(this->lookupOrDefault<scalar>("epsilon", 1.0e-10)),
    epsilonMin_(this->lookupOrDefault<scalar>("epsilonMin", epsilon_*0.1)),
    epsilonMax_(this->lookupOrDefault<scalar>("epsilonMax", epsilon_*100)),
    pthermo1_(rhoThermo::New(mesh_, "", "gasThermo")),
    pthermo2_(rhoThermo::New(mesh_, "", "liquidThermo")),
    thermo1_(pthermo1_()),
    thermo2_(pthermo2_()),
    pGasProps1_(gasProperties::New(thermo1_)),
    pGasProps2_(gasProperties::New(thermo2_)),
    gasProps1_(pGasProps1_()),
    gasProps2_(pGasProps2_()),
    p_(thermo1_.p()),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    T_(thermo1_.T()),
    w_
    (
        IOobject
        (
            "w",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimLength/dimTime
    ),
    conservative_
    (
        p_,
        U_,
        T_,
        w_,
        thermo1_,
        thermo2_
    )
{
    correctBoundaryCondition();
    correctThermo();
    correctConservative();
}


void Foam::compressibleMixture::correctSoundSpeeds()
{
    forAll(a_, celli)
    {
        a_[celli] = gasProps1_.c(p_[celli], T_[celli])*sqrt(1.0 - w_[celli]);
    }

    const volScalarField::Boundary& pBf = p_.boundaryField();
    const volScalarField::Boundary& TBf = T_.boundaryFieldRef();
    const volScalarField::Boundary& wBf = w_.boundaryFieldRef();
    volScalarField::Boundary& aBf = a_.boundaryFieldRef();
    
    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& pw = wBf[patchi];
        fvPatchScalarField& pa = aBf[patchi];

        forAll(pp, facei)
        {
            pa[facei] = gasProps1_.c(pp[facei], pT[facei])*sqrt(1.0 - pw[facei]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::compressibleMixture::correct()
{
    U_ = conservative_.rhoU()/conservative_.rho();
    w_ = conservative_.rhow()/conservative_.rho();

    //thermo1_.rho() = (1.0 - w_)*conservative_.rho();
    
    volScalarField e = conservative_.rhoE()/conservative_.rho() - 0.5*magSqr(U_);

    //Newton
    forAll(p_, celli)
    {
        const scalar tol = 1.e-6;
        const label maxIter = 100;
        
        scalar p = p_[celli];
        scalar T = T_[celli];

        const scalar pTol = p*tol;
        const scalar TTol = T*tol;
        
        const scalar w = w_[celli];
        const scalar rho1 = (1.0 - w)*conservative_.rho()[celli];
        //const scalar v = 1.0/conservative_.rho()[celli];

        label iter = 0;
        bool exitLoop = false;

        while(!exitLoop)
        {
            const scalar v1 = 1.0/gasProps1_.rho(p, T);
            const scalar v2 = 1.0/gasProps2_.rho(p, T);

            const scalar beta_T1 = gasProps1_.beta_T(p, T);
            const scalar beta_T2 = gasProps2_.beta_T(p, T);
            const scalar beta_p1 = gasProps1_.beta_p(p, T);
            const scalar beta_p2 = gasProps2_.beta_p(p, T);

            //const scalar f1 = (1.0 - w)*v1 + w*v2 - v;
            const scalar f1 = v1 - 1.0/rho1;
            const scalar f2 = (1.0 - w)*gasProps1_.Es(p, T) + w*gasProps2_.Es(p, T) - e[celli];
            
            //const scalar df1dp = -(1.0 - w)*v1*beta_T1 - w*v2*beta_T2;
            //const scalar df1dT =  (1.0 - w)*v1*beta_p1 + w*v2*beta_p2;

            const scalar df1dp = -v1*beta_T1;
            const scalar df1dT =  v1*beta_p1;
            const scalar df2dp = p*v2*beta_T2 - T*v2*beta_p2 + (1.0 - w)*p*v1*beta_T1 - (1.0 - w)*T*v1*beta_p2;
            const scalar df2dT = w*gasProps2_.Cp(p, T) - w*p*v2*beta_p2 + (1.0 - w)*gasProps1_.Cp(p, T) - (1.0 - w)*p*v1*beta_p1;

            const scalar dp = ( df2dT*f1 - df1dT*f2)/(df1dp*df2dT - df1dT*df2dp);
            const scalar dT = (-df2dp*f1 + df1dp*f2)/(df1dp*df2dT - df1dT*df2dp);

            p -= dp;
            T -= dT;

            exitLoop = (mag(dp) < pTol) && (mag(dT) < TTol);

            if (iter++ > maxIter)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter
                        << ", p : " << p
                        << ", T : " << T
                        << ", tol_p: " << mag(dp)
                        << ", tol_T: " << mag(dT)
                        << abort(FatalError);
            }
        }
        p_[celli] = p;
        T_[celli] = T;
    } //end Newton
}


void Foam::compressibleMixture::correctBoundaryCondition()
{
    p_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
    w_.correctBoundaryConditions();
}


void Foam::compressibleMixture::correctThermo()
{
    thermo1_.he() = thermo1_.he(p_, T_);
    thermo2_.he() = thermo2_.he(p_, T_);

    //TODO thermo correct update z (p, T) - nove "TRhoThermo.H"
    thermo1_.correct();
    thermo2_.correct();

    correctSoundSpeeds();
}

void Foam::compressibleMixture::correctConservative()
{
    conservative_.rho() = rho();
    conservative_.rhoU() = conservative_.rho()*U_;
    conservative_.rhoE() = (1.0 - w_)*thermo1_.he() + conservative_.rho()*w_*thermo2_.he() + 0.5*Foam::magSqr(U_);
    conservative_.rhow() = conservative_.rho()*w_;

    //nejspise nepotrebuji - neni blending - netreba korece
}

tmp<surfaceScalarField> Foam::compressibleMixture::amaxSf() const
{
    return mag(fvc::interpolate(U_) & mesh_.Sf()) + mesh_.magSf()*fvc::interpolate(a_);
}

// ************************************************************************* //