/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration    | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "godunovInterphase.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(godunovInterphase, 0);
addToRunTimeSelectionTable(riemannSolver, godunovInterphase, dict);

namespace
{

struct interfaceState
{
    scalar rho;
    scalar soundSpeed;
    scalar gamma;
    scalar pressureOffset;
    vector velocity;
    scalar normalVelocity;
};

struct waveData
{
    scalar velocityJump;
    scalar derivative;
};

interfaceState makeInterfaceState
(
    const scalar p,
    const scalar T,
    const vector& U,
    const vector& normal,
    const gasProperties& gas
)
{
    interfaceState state;

    state.rho = gas.rho(p, T);
    state.soundSpeed = gas.c(p, T);
    state.gamma = gas.Cp(p, T)/gas.Cv(p, T);
    state.pressureOffset =
        state.rho*sqr(state.soundSpeed)/state.gamma - p;
    state.velocity = U;
    state.normalVelocity = U & normal;

    if
    (
        state.rho <= VSMALL
     || state.soundSpeed <= VSMALL
     || state.gamma <= 1.0
    )
    {
        FatalErrorInFunction
            << "Godunov interphase states require positive density and sound "
            << "speed, and gamma > 1. Got rho=" << state.rho
            << ", c=" << state.soundSpeed << ", gamma=" << state.gamma
            << exit(FatalError);
    }

    return state;
}

waveData waveFunction
(
    const interfaceState& state,
    const scalar p,
    const scalar starPressure
)
{
    waveData result;
    const scalar shiftedPressure = p + state.pressureOffset;
    const scalar shiftedStarPressure =
        starPressure + state.pressureOffset;

    if (shiftedPressure <= VSMALL || shiftedStarPressure <= VSMALL)
    {
        FatalErrorInFunction
            << "Godunov stiffened-gas wave requires positive shifted "
            << "pressures. Got p+pInf=" << shiftedPressure
            << " and pStar+pInf=" << shiftedStarPressure
            << exit(FatalError);
    }

    if (starPressure > p)
    {
        const scalar ratio =
            (state.gamma - 1.0)/(state.gamma + 1.0);
        const scalar denominator =
            shiftedStarPressure + ratio*shiftedPressure;
        const scalar coefficient =
            sqrt(2.0/((state.gamma + 1.0)*state.rho*denominator));

        result.velocityJump =
            (starPressure - p)*coefficient;
        result.derivative =
            coefficient
           *(
                1.0
              - 0.5*(starPressure - p)/denominator
            );
    }
    else
    {
        const scalar exponent =
            (state.gamma - 1.0)/(2.0*state.gamma);
        const scalar pressureRatio =
            shiftedStarPressure/shiftedPressure;

        result.velocityJump =
            2.0*state.soundSpeed/(state.gamma - 1.0)
           *(pow(pressureRatio, exponent) - 1.0);
        result.derivative =
            state.soundSpeed/(state.gamma*shiftedPressure)
           *pow
            (
                pressureRatio,
                -(state.gamma + 1.0)/(2.0*state.gamma)
            );
    }

    return result;
}

scalar solveStarPressure
(
    const interfaceState& left,
    const interfaceState& right,
    const scalar pLeft,
    const scalar pRight,
    const scalar velocityTolerance,
    const label maxIterations,
    waveData& leftWave,
    waveData& rightWave
)
{
    const scalar pressureScale =
        max
        (
            max(mag(pLeft), mag(pRight)),
            max
            (
                max(mag(left.pressureOffset), mag(right.pressureOffset)),
                scalar(1)
            )
        );
    scalar lowerPressure =
        max(-left.pressureOffset, -right.pressureOffset)
      + SMALL*pressureScale;

    const scalar lowerResidual =
        left.normalVelocity
      - right.normalVelocity
      + waveFunction(left, pLeft, lowerPressure).velocityJump
      + waveFunction(right, pRight, lowerPressure).velocityJump;

    if (lowerResidual >= 0)
    {
        FatalErrorInFunction
            << "The exact stiffened-gas Riemann problem has no positive "
            << "star-pressure solution (vacuum state)."
            << exit(FatalError);
    }

    scalar upperPressure =
        max(max(pLeft, pRight), lowerPressure) + pressureScale;
    scalar upperResidual = 0;
    bool bracketed = false;

    for (label i = 0; i < maxIterations; ++i)
    {
        upperResidual =
            left.normalVelocity
          - right.normalVelocity
          + waveFunction(left, pLeft, upperPressure).velocityJump
          + waveFunction(right, pRight, upperPressure).velocityJump;

        if (upperResidual > 0)
        {
            bracketed = true;
            break;
        }

        upperPressure = max(2.0*upperPressure, upperPressure + pressureScale);
    }

    if (!bracketed)
    {
        FatalErrorInFunction
            << "Unable to bracket the exact Godunov star pressure after "
            << maxIterations << " iterations."
            << exit(FatalError);
    }

    scalar starPressure = min
    (
        max(0.5*(pLeft + pRight), lowerPressure),
        upperPressure
    );

    for (label i = 0; i < maxIterations; ++i)
    {
        leftWave = waveFunction(left, pLeft, starPressure);
        rightWave = waveFunction(right, pRight, starPressure);

        const scalar residual =
            left.normalVelocity
          - right.normalVelocity
          + leftWave.velocityJump
          + rightWave.velocityJump;

        if (mag(residual) <= velocityTolerance)
        {
            return starPressure;
        }

        if (residual < 0)
        {
            lowerPressure = starPressure;
        }
        else
        {
            upperPressure = starPressure;
        }

        const scalar derivative =
            leftWave.derivative + rightWave.derivative;
        const scalar newtonPressure =
            starPressure - residual/stabilise(derivative, VSMALL);

        starPressure =
            newtonPressure > lowerPressure
         && newtonPressure < upperPressure
          ? newtonPressure
          : 0.5*(lowerPressure + upperPressure);
    }

    FatalErrorInFunction
        << "Exact Godunov star pressure did not converge to velocity "
        << "tolerance " << velocityTolerance << " in " << maxIterations
        << " iterations."
        << exit(FatalError);

    return starPressure;
}

scalar starDensity
(
    const interfaceState& state,
    const scalar p,
    const scalar starPressure
)
{
    const scalar pressureRatio =
        (starPressure + state.pressureOffset)
       /(p + state.pressureOffset);

    if (starPressure > p)
    {
        const scalar h = (state.gamma + 1.0)/(state.gamma - 1.0);
        return state.rho*(h*pressureRatio + 1.0)/(h + pressureRatio);
    }

    return state.rho*pow(pressureRatio, 1.0/state.gamma);
}

} // End anonymous namespace


godunovInterphase::godunovInterphase()
:
    riemannSolver(),
    velocityTolerance_(1.0e-4),
    maxIterations_(100)
{}


godunovInterphase::godunovInterphase(const dictionary& dict)
:
    riemannSolver(dict),
    velocityTolerance_
    (
        dict.lookupOrDefault<scalar>("velocityTolerance", 1.0e-4)
    ),
    maxIterations_(dict.lookupOrDefault<label>("maxIterations", 100))
{
    if (velocityTolerance_ <= 0 || maxIterations_ <= 0)
    {
        FatalErrorInFunction
            << "velocityTolerance and maxIterations must be positive, got "
            << velocityTolerance_ << " and " << maxIterations_
            << exit(FatalError);
    }
}


void godunovInterphase::calculateFlux
(
    scalar& starPressure,
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
    starPressure = 0;
    alphaRhoFlux1Left = alphaRhoFlux1Right = 0;
    alphaRhoFlux2Left = alphaRhoFlux2Right = 0;
    alphaRhoUFlux1Left = alphaRhoUFlux1Right = vector::zero;
    alphaRhoUFlux2Left = alphaRhoUFlux2Right = vector::zero;
    alphaRhoEFlux1Left = alphaRhoEFlux1Right = 0;
    alphaRhoEFlux2Left = alphaRhoEFlux2Right = 0;

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
    const scalar alphaJump = alphaR - alphaL;
    const scalar interfaceLength = mag(alphaJump);

    if (interfaceLength <= VSMALL)
    {
        return;
    }

    const vector normal = Sf/magSf;
    const bool phase1OnLeft = alphaJump > 0;

    const interfaceState leftState = phase1OnLeft
      ? makeInterfaceState(pLeft, T1Left, U1Left, normal, gas1)
      : makeInterfaceState(pLeft, T2Left, U2Left, normal, gas2);
    const interfaceState rightState = phase1OnLeft
      ? makeInterfaceState(pRight, T2Right, U2Right, normal, gas2)
      : makeInterfaceState(pRight, T1Right, U1Right, normal, gas1);

    waveData leftWave{0, 0};
    waveData rightWave{0, 0};
    starPressure = solveStarPressure
    (
        leftState,
        rightState,
        pLeft,
        pRight,
        velocityTolerance_,
        maxIterations_,
        leftWave,
        rightWave
    );

    const scalar starNormalVelocity =
        0.5
       *(
            leftState.normalVelocity - leftWave.velocityJump
          + rightState.normalVelocity + rightWave.velocityJump
        );
    const bool leftPhaseDonates = starNormalVelocity > 0;
    const interfaceState& donorState =
        leftPhaseDonates ? leftState : rightState;
    const scalar donorPressure = leftPhaseDonates ? pLeft : pRight;
    const scalar donorStarDensity =
        starDensity(donorState, donorPressure, starPressure);
    const vector donorStarVelocity =
        donorState.velocity
      + (starNormalVelocity - donorState.normalVelocity)*normal;
    const scalar donorStarEnthalpy =
        donorState.gamma/(donorState.gamma - 1.0)
       *(starPressure + donorState.pressureOffset)
       /donorStarDensity
      + 0.5*magSqr(donorStarVelocity);

    const scalar interfaceMassFlux =
        interfaceLength*magSf*donorStarDensity*starNormalVelocity;
    const vector interfaceMomentumFlux =
        interfaceMassFlux*donorStarVelocity;
    const scalar interfaceEnergyFlux =
        interfaceMassFlux*donorStarEnthalpy;

    scalar& donorMassFluxLeft =
        leftPhaseDonates
      ? (phase1OnLeft ? alphaRhoFlux1Left : alphaRhoFlux2Left)
      : (phase1OnLeft ? alphaRhoFlux2Left : alphaRhoFlux1Left);
    scalar& donorMassFluxRight =
        leftPhaseDonates
      ? (phase1OnLeft ? alphaRhoFlux1Right : alphaRhoFlux2Right)
      : (phase1OnLeft ? alphaRhoFlux2Right : alphaRhoFlux1Right);
    vector& donorMomentumFluxLeft =
        leftPhaseDonates
      ? (phase1OnLeft ? alphaRhoUFlux1Left : alphaRhoUFlux2Left)
      : (phase1OnLeft ? alphaRhoUFlux2Left : alphaRhoUFlux1Left);
    vector& donorMomentumFluxRight =
        leftPhaseDonates
      ? (phase1OnLeft ? alphaRhoUFlux1Right : alphaRhoUFlux2Right)
      : (phase1OnLeft ? alphaRhoUFlux2Right : alphaRhoUFlux1Right);
    scalar& donorEnergyFluxLeft =
        leftPhaseDonates
      ? (phase1OnLeft ? alphaRhoEFlux1Left : alphaRhoEFlux2Left)
      : (phase1OnLeft ? alphaRhoEFlux2Left : alphaRhoEFlux1Left);
    scalar& donorEnergyFluxRight =
        leftPhaseDonates
      ? (phase1OnLeft ? alphaRhoEFlux1Right : alphaRhoEFlux2Right)
      : (phase1OnLeft ? alphaRhoEFlux2Right : alphaRhoEFlux1Right);

    donorMassFluxLeft += interfaceMassFlux;
    donorMassFluxRight += interfaceMassFlux;
    donorMomentumFluxLeft += interfaceMomentumFlux;
    donorMomentumFluxRight += interfaceMomentumFlux;
    donorEnergyFluxLeft += interfaceEnergyFlux;
    donorEnergyFluxRight += interfaceEnergyFlux;
}

} // End namespace Foam
