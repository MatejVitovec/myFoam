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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
/*Foam::TwoFluidFoam::twoFluid::twoFluid
(
    const volScalarField& p,
    const volScalarField& alpha,
    const volVectorField& U1,
    const volVectorField& U2,
    const volScalarField& T1,
    const volScalarField& T2,
    rhoThermo& thermo1,
    rhoThermo& thermo2
)
:
mesh_(mesh),
p_(p),
alpha_(alpha),
U1_(U1),
U2_(U2),
T1_(T1),
T2_(T2),
thermo1_(thermo1),
thermo2_(thermo2),
gasProps1_(gasProperties::New(thermo1)),
gasProps2_(gasProperties::New(thermo2)),
pInt_
(
    IOobject
    (
        "pInt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh_
),
twoFluidConservative_()
{
    
}*/

Foam::TwoFluidFoam::twoFluid::twoFluid(const fvMesh& mesh)
:
epsilon_(1.0e-7),
epsilonMin_(1.0e-1*epsilon_),
epsilonMax_(1.0e3*epsilon_),
mesh_(mesh),
pthermo1_(rhoThermo::New(mesh_, "1")),
pthermo2_(rhoThermo::New(mesh_, "2")),
thermo1_(pthermo1_()),
thermo2_(pthermo2_()),
pGasProps1_(gasProperties::New(thermo1_)),
pGasProps2_(gasProperties::New(thermo2_)),
gasProps1_(pGasProps1_()),
gasProps2_(pGasProps2_()),
p_(thermo1_.p()),
alpha_
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
    alpha_,
    U1_,
    U2_,
    T1_,
    T2_,
    thermo1_,
    thermo2_,
    pInt_
)
{
    //blendVanishingFluid();
    //correctBoundaryCondition();
    //correctThermo();
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
)
{
    const vector U1 = alphaRhoU1/alphaRho1;
    const vector U2 = alphaRhoU2/alphaRho2;

    const scalar U1magSqr = magSqr(U1);
    const scalar U2magSqr = magSqr(U2);

    //Info << "U1_old: " << U1Ref << " U1: " << U1 << endl;
    //Info << "U2_old: " << U2Ref << " U2: " << U2 << endl;
    //Info << "rho1: " << alphaRho1/alphaRef << endl;
    //Info << "rho2: " << alphaRho2/alphaRef << endl;

    // Estimate from old values
    scalar p = pRef;
    scalar T1 = T1Ref;
    scalar T2 = T2Ref;

    //Newton
    {
        const scalar tol = 1.e-6;
        const label maxIter = 100;

        const scalar TTol = T1*tol;
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

            exitLoop = (mag(dx[0]) < pTol) && (mag(dx[1]) < TTol) && (mag(dx[2]) < TTol);

            if (iter++ > maxIter)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter
                        //<< " when starting from p0:" << p0
                        << ", T1 : " << T1
                        << ", T2 : " << T2
                        << ", p : " << p
                        << ", T_base_tol: " << TTol
                        << ", tol_p: " << mag(dx[0])
                        << ", tol_T1: " << mag(dx[1])
                        << ", tol_T2: " << mag(dx[2])
                        << abort(FatalError);
            }
        }
    } //end Newton

    pRef = p;
    alphaRef = alphaRho1/gasProps1_.rho(p, T1);
    U1Ref = U1;
    U2Ref = U2;    
    T1Ref = T1;
    T2Ref = T2;    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TwoFluidFoam::twoFluid::correct()
{
    forAll(mesh_.cells(), celli)
    {
        primitiveFromConservative
        (
            p_[celli],
            alpha_[celli],
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
}

void Foam::TwoFluidFoam::twoFluid::blendVanishingFluid()
{
    forAll(mesh_.cells(), celli)
    {
        const scalar alpha = alpha_[celli];

        if ((1.0 - alpha) <= epsilonMin_)
        {
            alpha_[celli] = 1.0 - epsilonMin_;
            U2_[celli] = U1_[celli];
            T2_[celli] = T1_[celli];
        }
        else if ((1.0 - alpha) < epsilonMax_)
        {
            const scalar xi = ((1.0 - alpha) - epsilonMin_)/(epsilonMax_ - epsilonMin_);
            const scalar gFunc = -Foam::pow(xi, 2.0)*(2.0*xi - 3.0);
            
            U2_[celli] = gFunc*U2_[celli] + (1.0 - gFunc)*U1_[celli];
            T2_[celli] = gFunc*T2_[celli] + (1.0 - gFunc)*T1_[celli];
        }
        
        if (alpha <= epsilonMin_)
        {
            alpha_[celli] = epsilonMin_;
            U1_[celli] = U2_[celli];
            T1_[celli] = T2_[celli];
        }
        else if (alpha < epsilonMax_)
        {
            const double xi = (alpha - epsilonMin_)/(epsilonMax_ - epsilonMin_);
            const double gFunc = -std::pow(xi, 2.0)*(2.0*xi - 3.0);

            U1_[celli] = gFunc*U1_[celli] + (1.0 - gFunc)*U2_[celli];
            T1_[celli] = gFunc*T1_[celli] + (1.0 - gFunc)*T2_[celli];
        }
    }
}

void Foam::TwoFluidFoam::twoFluid::correctBoundaryCondition()
{
    p_.correctBoundaryConditions();
    alpha_.correctBoundaryConditions();
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
}

void Foam::TwoFluidFoam::twoFluid::correctInterfacialPressure()
{
    const scalar sigma = 2.0;
    const scalar epsilonP = 0.01;

    const auto& rho1 = thermo1_.rho();
    const auto& rho2 = thermo2_.rho();

    pInt_ = p_ - Foam::min(sigma*((alpha_*rho1*(1.0 - alpha_)*rho2)/
        (alpha_*rho2 + (1.0 - alpha_)*rho1))*Foam::magSqr(U2_ - U1_), epsilonP*p_);
}

void Foam::TwoFluidFoam::twoFluid::correctConservative()
{
    const auto& rho1 = thermo1_.rho();
    const auto& rho2 = thermo2_.rho();

    const auto E1 = thermo1_.he() + 0.5*Foam::magSqr(U1_);
    const auto E2 = thermo2_.he() + 0.5*Foam::magSqr(U2_);

    conservative_.alphaRho1() = alpha_*rho1;
    conservative_.alphaRho2() = (1.0 - alpha_)*rho2;
    conservative_.alphaRhoU1() = conservative_.alphaRho1()*U1_;
    conservative_.alphaRhoU2() = conservative_.alphaRho2()*U2_;
    conservative_.epsilon1() = alpha_*(rho1*E1 + pInt_);
    conservative_.epsilon2() = (1.0 - alpha_)*(rho2*E2 + pInt_);

    //nejspise nepotrebuji konzervativní pole staci pouze jako INTERNAL FIELD
    /*conservative_.alphaRho1().boundaryFieldRef() = alpha_.boundaryField()*rho1.boundaryField();
    conservative_.alphaRho2().boundaryFieldRef() = (1.0 - alpha_.boundaryField())*rho2.boundaryField();
    conservative_.alphaRhoU1().boundaryFieldRef() = alpha_.boundaryField()*rho1.boundaryField()*U1_.boundaryField();
    conservative_.alphaRhoU2().boundaryFieldRef() = (1.0 - alpha_.boundaryField())*rho2.boundaryField()*U2_.boundaryField();
    conservative_.epsilon1().boundaryFieldRef() = alpha_.boundaryField()*(rho1.boundaryField()*E1.boundaryField() + pInt_.boundaryField());
    conservative_.epsilon2().boundaryFieldRef() = (1.0 - alpha_.boundaryField())*(rho2.boundaryField()*E2.boundaryField() + pInt_.boundaryField());*/
}

// ************************************************************************* //