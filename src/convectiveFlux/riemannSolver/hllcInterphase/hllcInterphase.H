/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

Class
    hllcInterphase

Description
    HLLC flux contribution for an interface between two phases. Returns the
    HLLC star pressure for use by the composite solver's interfacial-pressure
    correction.

SourceFiles
    hllcInterphase.C

\*---------------------------------------------------------------------------*/

#ifndef hllcInterphase_H
#define hllcInterphase_H

#include "riemannSolver.H"

namespace Foam
{

class hllcInterphase
:
    public riemannSolver
{
    const scalar minimumStarPressureRatio_;

    hllcInterphase(const hllcInterphase&) = delete;
    void operator=(const hllcInterphase&) = delete;

    scalar einfeldtSoundSpeed
    (
        const scalar pLeft,    const scalar pRight,
        const vector ULeft,    const vector URight,
        const scalar TLeft,    const scalar TRight,
        const scalar rhoLeft,  const scalar rhoRight,
        const vector normal,
        const gasProperties& gas1,
        const gasProperties& gas2
    ) const;

    scalar roeSoundSpeed
    (
        const scalar pLeft,    const scalar pRight,
        const vector ULeft,    const vector URight,
        const scalar TLeft,    const scalar TRight,
        const scalar rhoLeft,  const scalar rhoRight,
        const vector normal,
        const gasProperties& gas1,
        const gasProperties& gas2
    ) const;

public:

    TypeName("hllcInterphase");

    hllcInterphase();
    explicit hllcInterphase(const dictionary& dict);

    virtual ~hllcInterphase() = default;

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
