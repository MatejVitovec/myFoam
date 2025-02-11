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
        rhoThermo::New(mesh)
    );
    thermo1_ = pThermo1();

    autoPtr<rhoThermo> pThermo2
    (
        rhoThermo::New(mesh)
    );
    thermo2_ = pThermo2();

    p_ = thermo1.p();
    T1_ = thermo1.T();
    T2_ = thermo2.T(); //TODO

    alpha_ = volScalarField
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
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




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::TwoFluidFoam::twoFluid::correct()
{
     
}

void Foam::TwoFluidFoam::twoFluid::blend()
{
    const scalar epsilonMin = 1.0e-5;
    const scalar epsilonMax = 1.0e-4;

    forAll(mesh.V(), celli)
    {
        const scalar alpha = alpha_[celli]

        if (1.0 - alpha <= epsilonMin)
        {
            alpha[celli] = 1.0 - epsilonMin;
            U2[celli] = U1[celli];
            T2[celli] = T1[celli];
        }
        else if (1.0 - alpha < epsilonMax)
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