/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "condensationMomentum.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationMomentum, 0);
addToRunTimeSelectionTable(condensationModel, condensationMomentum, params);


condensationMomentum::condensationMomentum
(
    volScalarField& alphaL,
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur,
    const dictionary& dict
)
:
    condensationModel(alphaL, alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, satur, dict),

    Q0_(
        IOobject
        (
            "Q0",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Q1_(
        IOobject
        (
            "Q1",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    Q2_(
        IOobject
        (
            "Q2",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    pNucleation_(nucleationModel::New(dict, gasThermo_, liquidThermo_, saturation_, gasProps_, liquidProps_)),
    nucleation_(pNucleation_()),

    pGrowth_(growthModel::New(dict, gasThermo_, liquidThermo_, saturation_)),
    growth_(pGrowth_())
{}


void condensationMomentum::constrainMoments()
{
    const volScalarField& J = nucleation_.J();

    forAll(alphaL_, i)
    {
        //TODO alphaMin/wmin
        if (alphaL_[i] <= 1e-29 && J[i] <= 0)
        {
            Q0_[i] = 0;
            Q1_[i] = 0;
            Q2_[i] = 0;
        }

        Q0_[i] = max(Q0_[i], 0);
        Q1_[i] = max(Q1_[i], 0);
        Q2_[i] = max(Q2_[i], 0);
    }
}


void condensationMomentum::correct()
{
    nucleation_.correct();
    growth_.correct(dropletRadius());

    const fvMesh& mesh = gasThermo_.p().mesh();
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;
    fields.add(Q0_);
    fields.add(Q1_);
    fields.add(Q2_);

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            alphaRhoPhi_,
            mesh.divScheme("div(alphaRhoPhi2,w_Q)")
        )
    );

    {
        fvScalarMatrix Q2Eqn
        (
            fvm::ddt(alpha_, rho_, Q2_)
            + mvConvection->fvmDiv(alphaRhoPhi_, Q2_)
            ==
            alpha_*sqr(nucleation_.rc())*nucleation_.J() + 2*alpha_*rho_*Q1_*growth_.rDot()
        );

        Q2Eqn.relax();
        Q2Eqn.solve();
    }

    {
        fvScalarMatrix Q1Eqn
        (
            fvm::ddt(alpha_, rho_, Q1_) 
            + mvConvection->fvmDiv(alphaRhoPhi_, Q1_)
            ==
            alpha_*nucleation_.rc()*nucleation_.J() + alpha_*rho_*Q0_*growth_.rDot()
        );

        Q1Eqn.relax();
        Q1Eqn.solve();
    }

    {
        fvScalarMatrix Q0Eqn
        (
            fvm::ddt(alpha_, rho_, Q0_) 
            + mvConvection->fvmDiv(alphaRhoPhi_, Q0_)
            ==
            alpha_*nucleation_.J()
        );

        Q0Eqn.relax();
        Q0Eqn.solve();
    }

    constrainMoments();

    //nucleationRateMassSource_ = 4.0/3.0*pi*pow3(rc)*J*rho_l;
    //growthRateMassSource_ = 4*pi*(alpha_*rho_*Q2_)*rDot*rho_l;
    //growthRateMassSource_ = 4*pi*alpha_*(alpha_*rho_l + (1.0 - alpha_)*rho_g)*Q2_*rDot*rho_l;
}


tmp<volScalarField> condensationMomentum::dropletRadius() const
{
    return sqrt(mag(Q2_)/(mag(Q0_ + dimensionedScalar("Q0SMALL", Q0_.dimensions(), SMALL))));
}


tmp<volScalarField> condensationMomentum::dropletDiameter() const
{
    return 2.0*max(dropletRadius(), growth_.rMin());
}


tmp<volScalarField> condensationMomentum::nucleationRateMassSource() const
{
    const scalar pi = constant::mathematical::pi;

    return 4/3*pi*pow3(nucleation_.rc())*nucleation_.J()*liquidThermo_.rho(); // TODO nevim jak s alphou
}


tmp<volScalarField> condensationMomentum::growthRateMassSource() const
{
    const scalar pi = constant::mathematical::pi;
    
    return 4*pi*alpha_*rho_*Q2_*growth_.rDot()*liquidThermo_.rho();
}

/*tmp<volScalarField> condensationMomentum::r32() const
{
    const scalar pi = constant::mathematical::pi;

    //TODO n

    return (3.0*alpha_)/(4.0*pi*n*rho_*Q2_)
}

tmp<volScalarField> condensationMomentum::r30() const
{
    const scalar pi = constant::mathematical::pi;

    //TODO n and y

    return pow((3.0*y)/(4*pi*liquidThermo_.rho()*n), 1.0/3.0);
}

tmp<volScalarField> condensationMomentum::r20() const
{
    return Q2_/(Q0_ + dimensionedScalar("Q0SMALL", Q0_.dimensions(), SMALL));
}

tmp<volScalarField> condensationMomentum::rG() const
{
    return r20()/(exp(0.5*sqr(log(sigmaG()))));
}

tmp<volScalarField> condensationMomentum::sigmaG() const
{
    return exp(sqrt(log(mag((Q0_*Q2_)/sqr(Q1_) - 1.0) + 1.0)));
}*/

}
}

// ************************************************************************* //
