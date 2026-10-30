/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration    |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "twoFluidComposite.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(twoFluidComposite, 0);
addToRunTimeSelectionTable(riemannSolver, twoFluidComposite, dict);


twoFluidComposite::twoFluidComposite(const dictionary& dict)
:
    riemannSolver(dict),
    interfaceTolerance_
    (
        dict.subDict("twoFluidCompositeCoeffs")
            .lookupOrDefault<scalar>("interfaceTolerance", 1.0e-5)
    ),
    interfacePressureCoefficient_
    (
        dict.subDict("twoFluidCompositeCoeffs")
            .lookupOrDefault<scalar>("interfacePressureCoefficient", 2.0)
    ),
    interfacePressureLimit_
    (
        dict.subDict("twoFluidCompositeCoeffs")
            .lookupOrDefault<scalar>("interfacePressureLimit", 0.01)
    ),
    phase1Solver_
    (
        riemannSolver::New
        (
            dict.subDict("twoFluidCompositeCoeffs").subDict("phase1Solver")
        )
    ),
    phase2Solver_
    (
        riemannSolver::New
        (
            dict.subDict("twoFluidCompositeCoeffs").subDict("phase2Solver")
        )
    ),
    interphaseSolver_
    (
        riemannSolver::New
        (
            dict.subDict("twoFluidCompositeCoeffs").subDict("interphaseSolver")
        )
    )
{
    if
    (
        interfaceTolerance_ < 0
     || interfacePressureCoefficient_ < 1
     || interfacePressureLimit_ < 0
    )
    {
        FatalErrorInFunction
            << "Invalid composite interphase coefficients: interfaceTolerance "
            << "and interfacePressureLimit must be non-negative, and "
            << "interfacePressureCoefficient must be at least one. Got "
            << interfaceTolerance_ << ", "
            << interfacePressureCoefficient_ << ", "
            << interfacePressureLimit_
            << exit(FatalError);
    }
}


vector twoFluidComposite::interphasePressureFlux
(
    const scalar alphaLeft,
    const scalar alphaRight,
    const scalar starPressure,
    const scalar rho1Interface,
    const scalar rho2Interface,
    const vector& U1Interface,
    const vector& U2Interface,
    const vector& Sf
) const
{
    const scalar interfaceLength = mag(alphaRight - alphaLeft);

    if (interfaceLength <= interfaceTolerance_)
    {
        return vector::zero;
    }

    const scalar alpha2Interface = 0.5*(alphaLeft + alphaRight);
    const scalar alpha1Interface = 1.0 - alpha2Interface;

    const scalar pressureDifference = min
    (
        interfacePressureCoefficient_
       *alpha2Interface*alpha1Interface*rho1Interface*rho2Interface
       /stabilise
        (
            alpha2Interface*rho1Interface
          + alpha1Interface*rho2Interface,
            VSMALL
        )
       *magSqr(U1Interface - U2Interface),
        interfacePressureLimit_*max(starPressure, VSMALL)
    );
    const scalar interfacePressure =
        max(starPressure - pressureDifference, VSMALL);

    return interfaceLength*interfacePressure*Sf;
}


void twoFluidComposite::calculateFlux
(
    scalar& alphaRhoFlux1Left,  scalar& alphaRhoFlux1Right,
    scalar& alphaRhoFlux2Left,  scalar& alphaRhoFlux2Right,
    vector& alphaRhoUFlux1Left, vector& alphaRhoUFlux1Right,
    vector& alphaRhoUFlux2Left, vector& alphaRhoUFlux2Right,
    scalar& alphaRhoEFlux1Left, scalar& alphaRhoEFlux1Right,
    scalar& alphaRhoEFlux2Left, scalar& alphaRhoEFlux2Right,
    const scalar& alphaLeft,    const scalar& alphaRight,
    const scalar& pLeft,        const scalar& pRight,
    const vector& U1Left,       const vector& U1Right,
    const vector& U2Left,       const vector& U2Right,
    const scalar& T1Left,       const scalar& T1Right,
    const scalar& T2Left,       const scalar& T2Right,
    const vector& Sf,
    const scalar& magSf,
    const gasProperties& gas1,
    const gasProperties& gas2
) const
{
    if (magSf <= VSMALL)
    {
        FatalErrorInFunction
            << "Face area magnitude must be positive, got " << magSf
            << exit(FatalError);
    }

    if
    (
        alphaLeft < -SMALL || alphaLeft > 1.0 + SMALL
     || alphaRight < -SMALL || alphaRight > 1.0 + SMALL
    )
    {
        FatalErrorInFunction
            << "Volume fractions must lie in [0, 1], got alphaLeft="
            << alphaLeft << " and alphaRight=" << alphaRight
            << exit(FatalError);
    }

    const scalar alphaL = min(max(alphaLeft, 0.0), 1.0);
    const scalar alphaR = min(max(alphaRight, 0.0), 1.0);

    if (mag(alphaR - alphaL) <= interfaceTolerance_)
    {
        const scalar commonSoundSpeed =
            0.25
           *(
                gas1.c(pLeft, T1Left)
              + gas1.c(pRight, T1Right)
              + gas2.c(pLeft, T2Left)
              + gas2.c(pRight, T2Right)
            );

        phase1Solver_->calculateFlux
        (
            alphaRhoFlux1Left,
            alphaRhoUFlux1Left,
            alphaRhoUFlux1Right,
            alphaRhoEFlux1Left,
            1.0 - alphaL,
            1.0 - alphaR,
            pLeft,
            pRight,
            U1Left,
            U1Right,
            gas1.rho(pLeft, T1Left),
            gas1.rho(pRight, T1Right),
            gas1.Es(pLeft, T1Left),
            gas1.Es(pRight, T1Right),
            commonSoundSpeed,
            commonSoundSpeed,
            Sf,
            magSf
        );
        alphaRhoFlux1Right = alphaRhoFlux1Left;
        alphaRhoEFlux1Right = alphaRhoEFlux1Left;

        phase2Solver_->calculateFlux
        (
            alphaRhoFlux2Left,
            alphaRhoUFlux2Left,
            alphaRhoUFlux2Right,
            alphaRhoEFlux2Left,
            alphaL,
            alphaR,
            pLeft,
            pRight,
            U2Left,
            U2Right,
            gas2.rho(pLeft, T2Left),
            gas2.rho(pRight, T2Right),
            gas2.Es(pLeft, T2Left),
            gas2.Es(pRight, T2Right),
            commonSoundSpeed,
            commonSoundSpeed,
            Sf,
            magSf
        );
        alphaRhoFlux2Right = alphaRhoFlux2Left;
        alphaRhoEFlux2Right = alphaRhoEFlux2Left;
        return;
    }

    const scalar alpha1Common = 1.0 - max(alphaL, alphaR);
    const scalar alpha2Common = min(alphaL, alphaR);

    scalar phase1MassFlux = 0;
    vector phase1MomentumFlux = vector::zero;
    scalar phase1EnergyFlux = 0;
    phase1Solver_->calculateFlux
    (
        phase1MassFlux,
        phase1MomentumFlux,
        phase1EnergyFlux,
        pLeft,
        pRight,
        U1Left,
        U1Right,
        T1Left,
        T1Right,
        Sf,
        magSf,
        gas1
    );

    scalar phase2MassFlux = 0;
    vector phase2MomentumFlux = vector::zero;
    scalar phase2EnergyFlux = 0;
    phase2Solver_->calculateFlux
    (
        phase2MassFlux,
        phase2MomentumFlux,
        phase2EnergyFlux,
        pLeft,
        pRight,
        U2Left,
        U2Right,
        T2Left,
        T2Right,
        Sf,
        magSf,
        gas2
    );

    alphaRhoFlux1Left = alphaRhoFlux1Right =
        alpha1Common*phase1MassFlux;
    alphaRhoUFlux1Left = alphaRhoUFlux1Right =
        alpha1Common*phase1MomentumFlux;
    alphaRhoEFlux1Left = alphaRhoEFlux1Right =
        alpha1Common*phase1EnergyFlux;

    alphaRhoFlux2Left = alphaRhoFlux2Right =
        alpha2Common*phase2MassFlux;
    alphaRhoUFlux2Left = alphaRhoUFlux2Right =
        alpha2Common*phase2MomentumFlux;
    alphaRhoEFlux2Left = alphaRhoEFlux2Right =
        alpha2Common*phase2EnergyFlux;

    scalar interfaceMassFlux1Left = 0;
    scalar interfaceMassFlux1Right = 0;
    scalar interfaceMassFlux2Left = 0;
    scalar interfaceMassFlux2Right = 0;
    vector interfaceMomentumFlux1Left = vector::zero;
    vector interfaceMomentumFlux1Right = vector::zero;
    vector interfaceMomentumFlux2Left = vector::zero;
    vector interfaceMomentumFlux2Right = vector::zero;
    scalar interfaceEnergyFlux1Left = 0;
    scalar interfaceEnergyFlux1Right = 0;
    scalar interfaceEnergyFlux2Left = 0;
    scalar interfaceEnergyFlux2Right = 0;
    scalar starPressure = 0;

    interphaseSolver_->calculateFlux
    (
        starPressure,
        interfaceMassFlux1Left,
        interfaceMassFlux1Right,
        interfaceMassFlux2Left,
        interfaceMassFlux2Right,
        interfaceMomentumFlux1Left,
        interfaceMomentumFlux1Right,
        interfaceMomentumFlux2Left,
        interfaceMomentumFlux2Right,
        interfaceEnergyFlux1Left,
        interfaceEnergyFlux1Right,
        interfaceEnergyFlux2Left,
        interfaceEnergyFlux2Right,
        alphaL,
        alphaR,
        pLeft,
        pRight,
        U1Left,
        U1Right,
        U2Left,
        U2Right,
        T1Left,
        T1Right,
        T2Left,
        T2Right,
        Sf,
        magSf,
        gas1,
        gas2
    );

    alphaRhoFlux1Left += interfaceMassFlux1Left;
    alphaRhoFlux1Right += interfaceMassFlux1Right;
    alphaRhoFlux2Left += interfaceMassFlux2Left;
    alphaRhoFlux2Right += interfaceMassFlux2Right;
    alphaRhoUFlux1Left += interfaceMomentumFlux1Left;
    alphaRhoUFlux1Right += interfaceMomentumFlux1Right;
    alphaRhoUFlux2Left += interfaceMomentumFlux2Left;
    alphaRhoUFlux2Right += interfaceMomentumFlux2Right;
    alphaRhoEFlux1Left += interfaceEnergyFlux1Left;
    alphaRhoEFlux1Right += interfaceEnergyFlux1Right;
    alphaRhoEFlux2Left += interfaceEnergyFlux2Left;
    alphaRhoEFlux2Right += interfaceEnergyFlux2Right;

    const bool phase1OnLeft = alphaR > alphaL;
    const scalar rho1Interface = phase1OnLeft
      ? gas1.rho(pLeft, T1Left)
      : gas1.rho(pRight, T1Right);
    const scalar rho2Interface = phase1OnLeft
      ? gas2.rho(pRight, T2Right)
      : gas2.rho(pLeft, T2Left);
    const vector U1Interface = phase1OnLeft ? U1Left : U1Right;
    const vector U2Interface = phase1OnLeft ? U2Right : U2Left;

    const vector pressureFlux = interphasePressureFlux
    (
        alphaL,
        alphaR,
        starPressure,
        rho1Interface,
        rho2Interface,
        U1Interface,
        U2Interface,
        Sf
    );

    if (phase1OnLeft)
    {
        alphaRhoUFlux1Left += pressureFlux;
        alphaRhoUFlux2Right += pressureFlux;
    }
    else
    {
        alphaRhoUFlux2Left += pressureFlux;
        alphaRhoUFlux1Right += pressureFlux;
    }
}

} // End namespace Foam
