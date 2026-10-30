/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#ifndef hllcGodunovInterphase_H
#define hllcGodunovInterphase_H

#include "riemannSolver.H"

namespace Foam
{

// HLLC estimates u* and p*; stiffened-gas Godunov relations provide the
// donor-side star state and interphase flux.
class hllcGodunovInterphase
:
    public riemannSolver
{
    const scalar minimumStarPressureRatio_;

    hllcGodunovInterphase(const hllcGodunovInterphase&) = delete;
    void operator=(const hllcGodunovInterphase&) = delete;

public:

    TypeName("hllcGodunovInterphase");

    hllcGodunovInterphase();
    explicit hllcGodunovInterphase(const dictionary& dict);

    virtual ~hllcGodunovInterphase() = default;

    virtual void calculateFlux
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
    ) const;
};

} // End namespace Foam

#endif
