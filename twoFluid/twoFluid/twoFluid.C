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
#include "twoFluid.H"
#include "boundMinMax.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
//TODO

Foam::TwoFluidFoam::twoFluid::twoFluid(const fvMesh& mesh)
:
IOdictionary
(
    IOobject
    (
        "twoFluidProperties",
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
pthermo1_(rhoThermo::New(mesh_, "1")),
pthermo2_(rhoThermo::New(mesh_, "2")),
thermo1_(pthermo1_()),
thermo2_(pthermo2_()),
pGasProps1_(gasProperties::New(thermo1_)),
pGasProps2_(gasProperties::New(thermo2_)),
gasProps1_(pGasProps1_()),
gasProps2_(pGasProps2_()),
p_(thermo1_.p()),
alpha2_
(
    IOobject
    (
        "alpha",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_
),
alpha1_
(
    IOobject
    (
        "alpha1",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1.0 - alpha2_
),
U1_
(
    IOobject
    (
        "U.1",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_
),
U2_
(
    IOobject
    (
        "U.2",
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh_
),
T1_(thermo1_.T()),
T2_(thermo2_.T()),
a1_
(
    IOobject
    (
        "a1",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimLength/dimTime
),
a2_
(
    IOobject
    (
        "a2",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimLength/dimTime
),
pInt_
(
    IOobject
    (
        "pInt",
        mesh_.time().timeName(),
        mesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    p_
),
conservative_
(
    p_,
    alpha2_,
    U1_,
    U2_,
    T1_,
    T2_,
    thermo1_,
    thermo2_,
    pInt_
)
{
    blendVanishingFluid();
    correctBoundaryCondition();
    correctThermo();
    correctInterfacialPressure();
    correctConservative();
}


void Foam::TwoFluidFoam::twoFluid::primitiveFromConservative
(
    scalar& pRef,
    scalar& alphaRef,
    vector& U1Ref,
    vector& U2Ref,
    scalar& T1Ref,
    scalar& T2Ref,
    const scalar alphaRho1,
    const scalar alphaRho2,
    const vector alphaRhoU1,
    const vector alphaRhoU2,
    const scalar epsilon1,
    const scalar epsilon2,
    const scalar pIntOld
) const
{
    const vector U1 = alphaRhoU1/alphaRho1;
    const vector U2 = alphaRhoU2/alphaRho2;

    const scalar U1magSqr = magSqr(U1);
    const scalar U2magSqr = magSqr(U2);

    // Estimate from old values
    scalar p = pRef;
    scalar T1 = T1Ref;
    scalar T2 = T2Ref;

    //Newton
    {
        const scalar tol = 1.e-6;
        const label maxIter = 100;

        const scalar TTol1 = T1*tol;
        const scalar TTol2 = T2*tol;
        const scalar pTol = p*tol;

        label iter = 0;
        bool exitLoop = false;

        while(!exitLoop)
        {
            const scalar v1 = 1.0/gasProps1_.rho(p, T1);
            const scalar v2 = 1.0/gasProps2_.rho(p, T2);

            const scalar beta_T1 = gasProps1_.beta_T(p, T1);
            const scalar beta_T2 = gasProps2_.beta_T(p, T2);
            const scalar beta_p1 = gasProps1_.beta_p(p, T1);
            const scalar beta_p2 = gasProps2_.beta_p(p, T2);

            const scalar de1dp = -v1*T1*beta_p1 + p*v1*beta_T1;
            const scalar de2dp = -v2*T2*beta_p2 + p*v2*beta_T2;
            const scalar de1dT = gasProps1_.Cp(p, T1) - p*v1*beta_p1;
            const scalar de2dT = gasProps2_.Cp(p, T2) - p*v2*beta_p2;

            vector f
            (
                alphaRho1*v1 + alphaRho2*v2 - 1.0,
                gasProps1_.Es(p, T1) + pIntOld*v1 + 0.5*U1magSqr - epsilon1/alphaRho1,
                gasProps2_.Es(p, T2) + pIntOld*v2 + 0.5*U2magSqr - epsilon2/alphaRho2
            );

            tensor J
            (
                alphaRho1*(-beta_T1*v1) + alphaRho2*(-beta_T2*v2),
                alphaRho1*beta_p1*v1,
                alphaRho2*beta_p2*v2,
                de1dp - beta_T1*v1*pIntOld,
                de1dT + beta_p1*v1*pIntOld,
                0.0,
                de2dp - beta_T2*v2*pIntOld,
                0.0,
                de2dT + beta_p2*v2*pIntOld
            );

            vector dx = inv(J) & f;

            p -= dx[0];
            T1 -= dx[1];
            T2 -= dx[2];

            exitLoop = (mag(dx[0]) < pTol) && (mag(dx[1]) < TTol1) && (mag(dx[2]) < TTol2);

            if (iter++ > maxIter)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter
                        //<< " when starting from p0:" << p0
                        << ", T1 : " << T1
                        << ", T2 : " << T2
                        << ", p : " << p
                        << ", T_1_tol: " << TTol1
                        << ", tol_p: " << mag(dx[0])
                        << ", tol_T1: " << mag(dx[1])
                        << ", tol_T2: " << mag(dx[2])
                        << abort(FatalError);
            }
        }
    } //end Newton

    pRef = p;
    //alphaRef = 1.0 - alphaRho1/gasProps1_.rho(p, T1);
    alphaRef = alphaRho2/gasProps2_.rho(p, T2);
    U1Ref = U1;
    U2Ref = U2;
    T1Ref = T1;
    T2Ref = T2;
}


void Foam::TwoFluidFoam::twoFluid::blendVanishingFluid
(
    scalar& alpha2,
    vector& U1,
    vector& U2,
    scalar& T1,
    scalar& T2
) const
{
    if ((1.0 - alpha2) <= epsilonMin_)
    {
        alpha2 = 1.0 - epsilonMin_;
        U1 = U2;
        T1 = T2;
    }
    else if ((1.0 - alpha2) < epsilonMax_)
    {
        const scalar xi = ((1.0 - alpha2) - epsilonMin_)/(epsilonMax_ - epsilonMin_);
        const scalar gFunc = -sqr(xi)*(2.0*xi - 3.0);
            
        U1 = gFunc*U1 + (1.0 - gFunc)*U2;
        T1 = gFunc*T1 + (1.0 - gFunc)*T2;
    }
        
    if (alpha2 <= epsilonMin_)
    {
        alpha2 = epsilonMin_;
        U2 = U1;
        T2 = T1;
    }
    else if (alpha2 < epsilonMax_)
    {
        const double xi = (alpha2 - epsilonMin_)/(epsilonMax_ - epsilonMin_);
        const double gFunc = -sqr(xi)*(2.0*xi - 3.0);

        U2 = gFunc*U2 + (1.0 - gFunc)*U1;
        T2 = gFunc*T2 + (1.0 - gFunc)*T1;
    }
}

void Foam::TwoFluidFoam::twoFluid::fluxFromConservative
(
    scalar& fluxAlphaRho1,
    scalar& fluxAlphaRho2,
    vector& fluxAlphaRhoU1,
    vector& fluxAlphaRhoU2,
    scalar& fluxAlphaRhoE1,
    scalar& fluxAlphaRhoE2,
    const scalar alphaRho1,
    const scalar alphaRho2,
    const vector alphaRhoU1,
    const vector alphaRhoU2,
    const scalar epsilon1,
    const scalar epsilon2,
    const scalar pIntOld,
    const vector Sf,
    const scalar p0,
    const scalar T10,
    const scalar T20
) const
{
    scalar p = p0;
    scalar alpha;
    vector U1;
    vector U2;
    scalar T1 = T10;
    scalar T2 = T20;

    primitiveFromConservative
    (
        p,
        alpha,
        U1,
        U2,
        T1,
        T2,
        alphaRho1,
        alphaRho2,
        alphaRhoU1,
        alphaRhoU2,
        epsilon1,
        epsilon2,
        pIntOld
    );

    //blending

    if ((1.0 - alpha) <= epsilonMin_)
    {
        alpha = 1.0 - epsilonMin_;
        U1 = U2;
        T1 = T2;

        const scalar alpha1 = 1.0 - alpha;
        const scalar rho1 = gasProps1_.rho(p, T1);
        const scalar E1 = gasProps1_.Es(p, T1) + 0.5*magSqr(U1);
        const scalar phi1 = U1 & Sf;
        fluxAlphaRho1 = alpha1*rho1*phi1;
        fluxAlphaRhoU1 = fluxAlphaRho1*U1 + alpha1*p*Sf;
        fluxAlphaRhoE1 = fluxAlphaRho1*E1;
    }
    else if ((1.0 - alpha) < epsilonMax_)
    {
        const scalar xi = ((1.0 - alpha) - epsilonMin_)/(epsilonMax_ - epsilonMin_);
        const scalar gFunc = -sqr(xi)*(2.0*xi - 3.0);
            
        U1 = gFunc*U1 + (1.0 - gFunc)*U2;
        T1 = gFunc*T1 + (1.0 - gFunc)*T2;

        const scalar alpha1 = 1.0 - alpha;
        const scalar rho1 = gasProps1_.rho(p, T1);
        const scalar E1 = gasProps1_.Es(p, T1) + 0.5*magSqr(U1);
        const scalar phi1 = U1 & Sf;
        fluxAlphaRho1 = alpha1*rho1*phi1;
        fluxAlphaRhoU1 = fluxAlphaRho1*U1 + alpha1*p*Sf;
        fluxAlphaRhoE1 = fluxAlphaRho1*E1;
    }
    else
    {
        const scalar alpha1 = 1.0 - alpha;
        const scalar phi1 = U1 & Sf;
        fluxAlphaRho1 = alphaRho1*phi1;
        fluxAlphaRhoU1 = fluxAlphaRho1*U1 + alpha1*p*Sf;
        fluxAlphaRhoE1 = ((epsilon1 - alpha1*pIntOld) + alpha1*p)*phi1;
    }
        
    if (alpha <= epsilonMin_)
    {
        alpha = epsilonMin_;
        U2 = U1;
        T2 = T1;

        const scalar rho2 = gasProps2_.rho(p, T2);
        const scalar E2 = gasProps2_.Es(p, T2) + 0.5*magSqr(U2);
        const scalar phi2 = U2 & Sf;
        fluxAlphaRho2 = alpha*rho2*phi2;
        fluxAlphaRhoU2 = fluxAlphaRho2*U2 + alpha*p*Sf;
        fluxAlphaRhoE2 = fluxAlphaRho2*E2;
    }
    else if (alpha < epsilonMax_)
    {
        const double xi = (alpha - epsilonMin_)/(epsilonMax_ - epsilonMin_);
        const double gFunc = -sqr(xi)*(2.0*xi - 3.0);

        U2 = gFunc*U2 + (1.0 - gFunc)*U1;
        T2 = gFunc*T2 + (1.0 - gFunc)*T1;

        const scalar rho2 = gasProps2_.rho(p, T2);
        const scalar E2 = gasProps2_.Es(p, T2) + 0.5*magSqr(U2);
        const scalar phi2 = U2 & Sf;
        fluxAlphaRho2 = alpha*rho2*phi2;
        fluxAlphaRhoU2 = fluxAlphaRho2*U2 + alpha*p*Sf;
        fluxAlphaRhoE2 = fluxAlphaRho2*E2;
    }
    else
    {
        const scalar phi2 = U2 & Sf;
        fluxAlphaRho2 = alphaRho2*phi2;
        fluxAlphaRhoU2 = fluxAlphaRho2*U2 + alpha*p*Sf;
        fluxAlphaRhoE2 = ((epsilon2 - alpha*pIntOld) + alpha*p)*phi2;
    }
}


void Foam::TwoFluidFoam::twoFluid::correctSoundSpeeds()
{
    forAll(a1_, celli)
    {
        a1_[celli] = gasProps1_.c(p_[celli], T1_[celli]);
        a2_[celli] = gasProps2_.c(p_[celli], T2_[celli]);
    }

    const volScalarField::Boundary& pBf = p_.boundaryField();
    const volScalarField::Boundary& T1Bf = T1_.boundaryFieldRef();
    const volScalarField::Boundary& T2Bf = T2_.boundaryFieldRef();
    volScalarField::Boundary& a1Bf = a1_.boundaryFieldRef();
    volScalarField::Boundary& a2Bf = a2_.boundaryFieldRef();
    
    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pT1 = T1Bf[patchi];
        const fvPatchScalarField& pT2 = T2Bf[patchi];
        fvPatchScalarField& pa1 = a1Bf[patchi];
        fvPatchScalarField& pa2 = a2Bf[patchi];

        forAll(pp, facei)
        {
            pa1[facei] = gasProps1_.c(pp[facei], pT1[facei]);
            pa2[facei] = gasProps2_.c(pp[facei], pT2[facei]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TwoFluidFoam::twoFluid::correct()
{
    forAll(mesh_.cells(), celli)
    {
        primitiveFromConservative
        (
            p_[celli],
            alpha2_[celli],
            U1_[celli],
            U2_[celli],
            T1_[celli],
            T2_[celli],
            conservative_.alphaRho1()[celli],
            conservative_.alphaRho2()[celli],
            conservative_.alphaRhoU1()[celli],
            conservative_.alphaRhoU2()[celli],
            conservative_.epsilon1()[celli],
            conservative_.epsilon2()[celli],
            pInt_[celli]
        );
    }

    alpha1_ = 1.0 - alpha2_;
}

void Foam::TwoFluidFoam::twoFluid::bound()
{
    const dimensionedScalar pMin = dimensionedScalar("pMin", dimPressure, 1e3);
    const dimensionedScalar pMax = dimensionedScalar("pMax", dimPressure, 1e12);
    const dimensionedScalar T1Min = dimensionedScalar("T1Min", dimTemperature, 50);
    const dimensionedScalar T1Max = dimensionedScalar("T1Max", dimTemperature, 1000);
    const dimensionedScalar T2Min = dimensionedScalar("T2Min", dimTemperature, 50);
    const dimensionedScalar T2Max = dimensionedScalar("T2Max", dimTemperature, 1000);

    boundMinMax(p_, pMin, pMax);
    boundMinMax(T1_, T1Min, T1Max);
    boundMinMax(T2_, T2Min, T2Max);
}

void Foam::TwoFluidFoam::twoFluid::blendVanishingFluid()
{
    forAll(mesh_.cells(), celli)
    {
        const scalar alpha2 = alpha2_[celli];

        if ((1.0 - alpha2) <= epsilonMin_)
        {
            alpha2_[celli] = 1.0 - epsilonMin_;
            U1_[celli] = U2_[celli];
            T1_[celli] = T2_[celli];
        }
        else if ((1.0 - alpha2) < epsilonMax_)
        {
            const scalar xi = ((1.0 - alpha2) - epsilonMin_)/(epsilonMax_ - epsilonMin_);
            const scalar gFunc = -sqr(xi)*(2.0*xi - 3.0);
            
            U1_[celli] = gFunc*U1_[celli] + (1.0 - gFunc)*U2_[celli];
            T1_[celli] = gFunc*T1_[celli] + (1.0 - gFunc)*T2_[celli];
        }
        
        if (alpha2 <= epsilonMin_)
        {
            alpha2_[celli] = epsilonMin_;
            U2_[celli] = U1_[celli];
            T2_[celli] = T1_[celli];
        }
        else if (alpha2 < epsilonMax_)
        {
            const double xi = (alpha2 - epsilonMin_)/(epsilonMax_ - epsilonMin_);
            const double gFunc = -sqr(xi)*(2.0*xi - 3.0);

            U2_[celli] = gFunc*U2_[celli] + (1.0 - gFunc)*U1_[celli];
            T2_[celli] = gFunc*T2_[celli] + (1.0 - gFunc)*T1_[celli];
        }
    }
    
    alpha1_ = 1.0 - alpha2_;
}

void Foam::TwoFluidFoam::twoFluid::correctBoundaryCondition()
{
    p_.correctBoundaryConditions();
    alpha2_.correctBoundaryConditions();
    alpha1_ = 1.0 - alpha2_;
    //alpha1_.boundaryFieldRef() = 1.0 - alpha2_.boundaryField();
    U1_.correctBoundaryConditions();
    U2_.correctBoundaryConditions();
    T1_.correctBoundaryConditions();
    T2_.correctBoundaryConditions();
}


void Foam::TwoFluidFoam::twoFluid::correctThermo()
{
    thermo1_.he() = thermo1_.he(p_, T1_);
    thermo2_.he() = thermo2_.he(p_, T2_);

    //TODO thermo correct update z (p, T) - nove "TRhoThermo.H"
    thermo1_.correct();
    thermo2_.correct();

    correctSoundSpeeds();
}

void Foam::TwoFluidFoam::twoFluid::correctInterfacialPressure()
{
    const scalar sigma = 2.0;
    const scalar epsilonP = 0.01;

    const auto& rho1 = thermo1_.rho();
    const auto& rho2 = thermo2_.rho();

    pInt_ = p_ - min(sigma*((alpha1_*rho1*alpha2_*rho2)/
        (alpha1_*rho2 + alpha2_*rho1))*magSqr(U2_ - U1_), epsilonP*p_);

    //pInt_ = p_ - min(sigma*alpha2_*rho1*magSqr(U2_ - U1_), epsilonP*p_);
}

void Foam::TwoFluidFoam::twoFluid::correctConservative()
{
    const auto& rho1 = thermo1_.rho();
    const auto& rho2 = thermo2_.rho();

    const auto E1 = thermo1_.he() + 0.5*Foam::magSqr(U1_);
    const auto E2 = thermo2_.he() + 0.5*Foam::magSqr(U2_);

    conservative_.alphaRho1() = alpha1_*rho1;
    conservative_.alphaRho2() = alpha2_*rho2;
    conservative_.alphaRhoU1() = conservative_.alphaRho1()*U1_;
    conservative_.alphaRhoU2() = conservative_.alphaRho2()*U2_;
    conservative_.epsilon1() = alpha1_*(rho1*E1 + pInt_);
    conservative_.epsilon2() = alpha2_*(rho2*E2 + pInt_);

    //nejspise nepotrebuji konzervativní pole staci pouze jako INTERNAL FIELD
    /*conservative_.alphaRho1().boundaryFieldRef() = alpha1_.boundaryField()*rho1.boundaryField();
    conservative_.alphaRho2().boundaryFieldRef() = (1.0 - alpha1_.boundaryField())*rho2.boundaryField();
    conservative_.alphaRhoU1().boundaryFieldRef() = alpha1_.boundaryField()*rho1.boundaryField()*U1_.boundaryField();
    conservative_.alphaRhoU2().boundaryFieldRef() = (1.0 - alpha1_.boundaryField())*rho2.boundaryField()*U2_.boundaryField();
    conservative_.epsilon1().boundaryFieldRef() = alpha1_.boundaryField()*(rho1.boundaryField()*E1.boundaryField() + pInt_.boundaryField());
    conservative_.epsilon2().boundaryFieldRef() = (1.0 - alpha1_.boundaryField())*(rho2.boundaryField()*E2.boundaryField() + pInt_.boundaryField());*/
}

tmp<surfaceScalarField> Foam::TwoFluidFoam::twoFluid::amaxSf() const
{
    return max(mag(fvc::interpolate(U1_) & mesh_.Sf()), mag(fvc::interpolate(U2_) & mesh_.Sf()))
         + max(mesh_.magSf()*fvc::interpolate(a1_), mesh_.magSf()*fvc::interpolate(a2_));
}

// ************************************************************************* //