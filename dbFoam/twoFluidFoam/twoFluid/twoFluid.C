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
Foam::TwoFluidFoam::twoFluid::twoFluid
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
p_(p);,
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
    
}

Foam::TwoFluidFoam::twoFluid::twoFluid(const fvMesh& mesh)
:
mesh_(mesh)
{
    autoPtr<rhoThermo> pThermo1
    (
        rhoThermo::New(mesh, "fluid1")
    );
    thermo1_ = pThermo1();

    autoPtr<rhoThermo> pThermo2
    (
        rhoThermo::New(mesh, "fluid2")
    );
    thermo2_ = pThermo2();

    p_ = thermo1.p();
    T1_ = thermo1.T();
    T2_ = thermo2.T();

    alpha_ = volScalarField
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );

    U1 = volVectorField
    (
        IOobject
        (
            "U1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );

    U2 = volVectorField
    (
        IOobject
        (
            "U2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    );

    pInt_ = volScalarField
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
    );

    twoFluidConservative_ = twoFluidConservative(p, alpha, U1, U2, T1, T2, thermo1, thermo2);
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
    gasProperties& gas1 = this->gas1();
    gasProperties& gas2 = this->gas2();

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

        const scalar TTol = T*tol;
        const scalar pTol = p*tol;

        label iter = 0;
        bool exitLoop = false;

        while(!exitLoop)
        {
            const scalar v1 = 1.0/gas1.rho(p, T1);
            const scalar v2 = 1.0/gas2.rho(p, T2);

            const scalar beta1_T = gas1.beta_T(p, T1);
            const scalar beta2_T = gas2.beta_T(p, T2);
            const scalar beta1_p = gas1.beta_p(p, T1);
            const scalar beta2_p = gas2.beta_p(p, T2);

            const scalar de1dp = -v1*T1*beta1_p + p*v1*beta1_T;
            const scalar de2dp = -v2*T2*beta2_p + p*v2*beta2_T;
            const scalar de1dT = gas1.Cp(p, T1) - p*v1*beta1_p;
            const scalar de2dT = gas2.Cp(p, T2) - p*v2*beta2_p;

            vector f
            (
                alphaRho1*v1 + alphaRho2*v2 - 1.0,
                alphaRho1*(gas1.Es(p, T1) + 0.5*U1magSqr + pIntOld*v1) - epsilon1,
                alphaRho2*(gas1.Es(p, T1) + 0.5*U2magSqr + pIntOld*v2) - epsilon2
            );

            tensor J
            (
                alphaRho1*(-beta1_T*v1) + alphaRho2*(-beta2_T*v2),
                alphaRho1*beta1_p*v1;,
                alphaRho2*beta2_p*v2;,
                alphaRho1*(1.0 + pIntOld)*de1dp,
                alphaRho1*(1.0 + pIntOld)*de1dT;,
                0.0,
                alphaRho2*(1.0 + pIntOld)*de2dp;,
                0.0,
                alphaRho2*(1.0 + pIntOld)*de2dT;
            );

            vector dx = inv(J) & f;

            p -= dx[0];
            T1 -= dx[1];
            T2 -= dx[2];

            exitLoop = (mag(dx[0]) < pTol) && (mag(dx[1]) < TTol) && (mag(dx[2]) < TTol)

            if (iter++ > maxIter)
            {
                FatalErrorInFunction
                    << "Maximum number of iterations exceeded: " << maxIter
                        //<< " when starting from p0:" << p0
                        << " T1  : " << T1
                        << " T2  : " << T2
                        << " p  : " << p
                        << " tol: " << tol
                        << abort(FatalError);
            }
        }
    } //end Newton

    pRef = p;
    alphaRef = alphaRho1/gas1.rho(p, T1);
    U1Ref = U1;
    U2Ref = U2;    
    T1Ref = T1;
    T2Ref = T2;    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TwoFluidFoam::twoFluid::correct()
{
    forAll(mesh.V(), celli)
    {
        primitiveFromConservative
        (
            p_[celli],
            alpha_[celli],
            U1_[celli],
            U2_[celli],
            T1_[celli],
            T2_[celli],
            alphaRho1_[celli],
            alphaRho2_[celli],
            alphaRhoU1_[celli],
            alphaRhoU2_[celli],
            epsilon1_[celli],
            epsilon2_[celli],
            pInt_[celli]
        )
    }
}

void Foam::TwoFluidFoam::twoFluid::blendVanishingFluid()
{
    const scalar epsilonMin = 1.0e-5;
    const scalar epsilonMax = 1.0e-4;

    forAll(mesh.V(), celli)
    {
        const scalar alpha = alpha_[celli]

        if ((1.0 - alpha) <= epsilonMin)
        {
            alpha[celli] = 1.0 - epsilonMin;
            U2[celli] = U1[celli];
            T2[celli] = T1[celli];
        }
        else if ((1.0 - alpha) < epsilonMax)
        {
            const scalar xi = ((1.0 - alpha_) - epsilonMin)/(epsilonMax - epsilonMin);
            const scalar gFunc = -Foam::pow(xi, 2.0)*(2.0*xi - 3.0);
            
            U2[celli] = gFunc*U2[celli] + (1.0 - gFunc)*U1[celli];
            T2[celli] = gFunc*T2[celli] + (1.0 - gFunc)*T1[celli];
        }
        
        if (alpha <= epsilonMin)
        {
            alpha[celli] = epsilonMin;
            U1[celli] = U2[celli];
            T1[celli] = T2[celli];
        }
        else if (alpha < epsilonMax)
        {
            const double xi = (alpha_[celli] - epsilonMin)/(epsilonMax - epsilonMin);
            const double gFunc = -std::pow(xi, 2.0)*(2.0*xi - 3.0);

            U1[celli] = gFunc*U1[celli] + (1.0 - gFunc)*U2[celli];
            T1[celli] = gFunc*T1[celli] + (1.0 - gFunc)*T2[celli];
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
    thermo1.he() = thermo1.he(p, T1);
    thermo2.he() = thermo2.he(p, T2);

    //TODO thermo correct update z (p, T) - nove "TRhoThermo.H"
    thermo1.correct();
    thermo2.correct();
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

    const auto& E1 = thermo1_.he() + 0.5*Foam::magSqr(U1);
    const auto& E2 = thermo2_.he() + 0.5*Foam::magSqr(U2)

    conservative_.alphaRho1() = alpha_*rho1;
    conservative_.alphaRho2() = (1.0 - alpha_)*rho2;
    conservative_.alphaRhoU1() = conservative_.alphaRho1()*U1_;
    conservative_.alphaRhoU2() = conservative_.alphaRho2()*U2_;
    conservative_.epsilon1() = alpha_*(rho1*E1 + pInt_);
    conservative_.epsilon2() = (1.0 - alpha_)*(rho2*E2 + pInt_);

    //nejspise nepotrebuji konzervativní pole staci pouze jako INTERNAL FIELD
    conservative_.alphaRho1().boundaryFieldRef() = alpha_.boundaryField()*rho1.boundaryField();
    conservative_.alphaRho2().boundaryFieldRef() = (1.0 - alpha_.boundaryField())*rho2.boundaryField();
    conservative_.alphaRhoU1().boundaryFieldRef() = alpha_.boundaryField()*rho1.boundaryField()*U1_.boundaryField();
    conservative_.alphaRhoU2().boundaryFieldRef() = (1.0 - alpha_.boundaryField())*rho2.boundaryField()*U2_.boundaryField();
    conservative_.epsilon1().boundaryFieldRef() = alpha_.boundaryField()*(rho1.boundaryField()*E1.boundaryField() + pInt_.boundaryField());
    conservative_.epsilon2().boundaryFieldRef() = (1.0 - alpha_.boundaryField())*(rho2.boundaryField()*E2.boundaryField() + pInt_.boundaryField());
}

// ************************************************************************* //