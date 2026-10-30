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
    twoFluidComposite

Description
    Composite two-fluid solver. It owns one single-fluid solver for each
    phase, a configurable solver for interphase flux, and the interfacial
    pressure correction.

    Configure the nested solvers in fvSchemes with, for example:

        riemannSolver twoFluidComposite;
        twoFluidCompositeCoeffs
        {
            phase1Solver { riemannSolver ausmPlusUp; }
            phase2Solver { riemannSolver ausmPlusUp; }
            interphaseSolver
            {
                riemannSolver godunovInterphase;
                velocityTolerance 1e-4;
                maxIterations 100;
            }
        }

    The phase solvers must implement the single-fluid calculateFlux overload
    and the volume-fraction-weighted overload used when alpha is continuous
    across a face. Advection-splitting solvers provide both overloads.

SourceFiles
    twoFluidComposite.C

\*---------------------------------------------------------------------------*/

#ifndef twoFluidComposite_H
#define twoFluidComposite_H

#include "riemannSolver.H"

namespace Foam
{

class twoFluidComposite
:
    public riemannSolver
{
    const scalar interfaceTolerance_;
    const scalar interfacePressureCoefficient_;
    const scalar interfacePressureLimit_;
    autoPtr<riemannSolver> phase1Solver_;
    autoPtr<riemannSolver> phase2Solver_;
    autoPtr<riemannSolver> interphaseSolver_;

    vector interphasePressureFlux
    (
        const scalar alphaLeft,
        const scalar alphaRight,
        const scalar starPressure,
        const scalar rho1Interface,
        const scalar rho2Interface,
        const vector& U1Interface,
        const vector& U2Interface,
        const vector& Sf
    ) const;

    twoFluidComposite(const twoFluidComposite&) = delete;
    void operator=(const twoFluidComposite&) = delete;

public:

    TypeName("twoFluidComposite");

    explicit twoFluidComposite(const dictionary& dict);

    virtual ~twoFluidComposite() = default;

    virtual void calculateFlux
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
    ) const;
};

} // End namespace Foam

#endif
