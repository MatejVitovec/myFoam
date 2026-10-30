/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /    O peration    |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "hllcInterphase.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

defineTypeNameAndDebug(hllcInterphase, 0);
addToRunTimeSelectionTable(riemannSolver, hllcInterphase, dict);

namespace
{

struct interfaceState
{
    scalar rho;
    scalar energy;
    scalar soundSpeed;
    scalar gamma;
    scalar pressureOffset;
    vector velocity;
    scalar normalVelocity;
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
    state.energy = gas.Es(p, T);
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
            << "HLLC interphase states require positive density and sound "
            << "speed, and gamma > 1.  Got rho=" << state.rho
            << ", c=" << state.soundSpeed << ", gamma=" << state.gamma
            << exit(FatalError);
    }

    return state;
}

void hllcFlux
(
    scalar& massFlux,
    vector& momentumFlux,
    scalar& energyFlux,
    const interfaceState& state,
    const scalar p,
    const scalar waveSpeed,
    const scalar contactSpeed,
    const vector& normal
)
{
    const scalar starDensity =
        state.rho*(waveSpeed - state.normalVelocity)
       /stabilise(waveSpeed - contactSpeed, VSMALL);

    const vector starVelocity =
        state.velocity
      + (contactSpeed - state.normalVelocity)*normal;

    const scalar specificTotalEnergy =
        state.energy + 0.5*magSqr(state.velocity);

    const scalar starSpecificTotalEnergy =
        specificTotalEnergy
      + (contactSpeed - state.normalVelocity)
       *(
            contactSpeed
          + p
           /stabilise
            (
                state.rho*(waveSpeed - state.normalVelocity),
                VSMALL
            )
        );

    const scalar stateEnergy = state.rho*specificTotalEnergy;
    const scalar starEnergy = starDensity*starSpecificTotalEnergy;

    massFlux =
        state.rho*state.normalVelocity
      + waveSpeed*(starDensity - state.rho);

    // The interfacial pressure contribution is added separately below.
    momentumFlux = massFlux*starVelocity;

    energyFlux =
        (stateEnergy + p)*state.normalVelocity
      + waveSpeed*(starEnergy - stateEnergy);
}

} // End anonymous namespace


hllcInterphase::hllcInterphase()
:
    riemannSolver(),
    minimumStarPressureRatio_(1.0e-6)
{}


hllcInterphase::hllcInterphase(const dictionary& dict)
:
    riemannSolver(dict),
    minimumStarPressureRatio_
    (
        dict.lookupOrDefault<scalar>("minimumStarPressureRatio", 1.0e-6)
    )
{
    if (minimumStarPressureRatio_ <= 0)
    {
        FatalErrorInFunction
            << "minimumStarPressureRatio must be positive, got "
            << minimumStarPressureRatio_
            << exit(FatalError);
    }
}


void hllcInterphase::calculateFlux
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

    const vector normal = Sf/magSf;
    const bool phase1OnLeft = alphaJump > 0;

    const interfaceState leftState = phase1OnLeft
      ? makeInterfaceState(pLeft, T1Left, U1Left, normal, gas1)
      : makeInterfaceState(pLeft, T2Left, U2Left, normal, gas2);
    const interfaceState rightState = phase1OnLeft
      ? makeInterfaceState(pRight, T2Right, U2Right, normal, gas2)
      : makeInterfaceState(pRight, T1Right, U1Right, normal, gas1);
    const scalar pInfLeft = leftState.pressureOffset;
    const scalar pInfRight = rightState.pressureOffset;

    const scalar sqrtRhoLeft = sqrt(leftState.rho);
    const scalar sqrtRhoRight = sqrt(rightState.rho);
    const scalar sqrtRhoSum = sqrtRhoLeft + sqrtRhoRight;
    const scalar inverseSqrtRhoSum = 1.0/sqrtRhoSum;

    const scalar gammaRoe =
        (
            sqrtRhoLeft*leftState.gamma
          + sqrtRhoRight*rightState.gamma
        )*inverseSqrtRhoSum;

    const vector velocityRoe =
        (
            sqrtRhoLeft*leftState.velocity
          + sqrtRhoRight*rightState.velocity
        )*inverseSqrtRhoSum;

    const scalar pressureDensityRoe =
        (gammaRoe - 1.0)/gammaRoe
       *(
            (
                sqrtRhoLeft*(pLeft + pInfLeft)/leftState.rho
              + sqrtRhoRight*(pRight + pInfRight)/rightState.rho
            )*inverseSqrtRhoSum
          + 0.5*sqrtRhoLeft*sqrtRhoRight*sqr(inverseSqrtRhoSum)
           *magSqr(leftState.velocity - rightState.velocity)
        );

    const scalar soundSpeedRoe =
        sqrt(max(gammaRoe*pressureDensityRoe, VSMALL));

    const scalar bLeft = min
    (
        min
        (
            leftState.normalVelocity - leftState.soundSpeed,
            (velocityRoe & normal) - soundSpeedRoe
        ),
        0.0
    );
    const scalar bRight = max
    (
        max
        (
            rightState.normalVelocity + rightState.soundSpeed,
            (velocityRoe & normal) + soundSpeedRoe
        ),
        0.0
    );

    const scalar contactSpeed =
        (
            pRight - pLeft
          + leftState.rho*leftState.normalVelocity
           *(bLeft - leftState.normalVelocity)
          - rightState.rho*rightState.normalVelocity
           *(bRight - rightState.normalVelocity)
        )
       /stabilise
        (
            leftState.rho*(bLeft - leftState.normalVelocity)
          - rightState.rho*(bRight - rightState.normalVelocity),
            VSMALL
        );

    const scalar starPressureLeft =
        pLeft
      + leftState.rho*(bLeft - leftState.normalVelocity)
       *(contactSpeed - leftState.normalVelocity);
    const scalar starPressureRight =
        pRight
      + rightState.rho*(bRight - rightState.normalVelocity)
       *(contactSpeed - rightState.normalVelocity);
    starPressure = max
    (
        0.5*(starPressureLeft + starPressureRight),
        minimumStarPressureRatio_*max(min(pLeft, pRight), VSMALL)
    );

    scalar massFlux = 0;
    vector momentumFlux = vector::zero;
    scalar energyFlux = 0;
    const bool leftPhaseDonates = contactSpeed > 0;

    if (leftPhaseDonates)
    {
        hllcFlux
        (
            massFlux,
            momentumFlux,
            energyFlux,
            leftState,
            pLeft,
            bLeft,
            contactSpeed,
            normal
        );
    }
    else
    {
        hllcFlux
        (
            massFlux,
            momentumFlux,
            energyFlux,
            rightState,
            pRight,
            bRight,
            contactSpeed,
            normal
        );
    }

    const scalar scaledInterfaceLength = interfaceLength*magSf;
    massFlux *= scaledInterfaceLength;
    momentumFlux *= scaledInterfaceLength;
    energyFlux *= scaledInterfaceLength;

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

    donorMassFluxLeft += massFlux;
    donorMassFluxRight += massFlux;
    donorMomentumFluxLeft += momentumFlux;
    donorMomentumFluxRight += momentumFlux;
    donorEnergyFluxLeft += energyFlux;
    donorEnergyFluxRight += energyFlux;

}

} // End namespace Foam
