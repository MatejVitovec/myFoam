/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "SchillerNaumann.H"
#include "twoFluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
    defineTypeNameAndDebug(SchillerNaumann, 0);
    addToRunTimeSelectionTable(dragModel, SchillerNaumann, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::SchillerNaumann::SchillerNaumann
(
    const dictionary& dict,
    const twoFluid& fluid,
    const volScalarField& dropletDiameter,
    const bool registerObject
)
:
    dragModel(dict, fluid, dropletDiameter, registerObject),
    residualRe_("residualRe", dimless, dict),
    Re_(
        IOobject
        (
            "Re",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    nu_(
        IOobject
        (
            "nu",
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.thermo1().nu()
    )
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::SchillerNaumann::~SchillerNaumann()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::TwoFluidFoam::SchillerNaumann::dCdRedx(const scalar Re, const scalar dRedx) const
{
    //neg(Re_ - 1000)*24.0*(1.0 + 0.15*pow(Re_, 0.687)) + pos0(Re_ - 1000)*0.44*max(Re_, residualRe_);
    return neg(Re - 1000)*((6183.0*dRedx)/(2500*pow(Re, 0.313))) + pos0(Re - 1000)*0.44*dRedx;
}


Foam::vector Foam::TwoFluidFoam::SchillerNaumann::dCdRedx(const scalar Re, const vector dRedx) const
{
    //return (6183.0*dRedx)/(2500*pow(Re, 0.313) + VSMALL);
    return neg(Re - 1000)*((6183.0*dRedx)/(2500*pow(Re, 0.313) + VSMALL)) + pos0(Re - 1000)*0.44*dRedx;
}


/*Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::SchillerNaumann::Cd() const
{
    volScalarField Re = max(mag(fluid_.U1() - fluid_.U2())*d_/fluid_.thermo1().nu(), 1e-6);

    //return neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687)) + pos0(Re - 1000)*0.44*max(Re, residualRe_);
    return neg(Re - 1000)*24.0*(1.0 + 0.15*pow(Re, 0.687))/Re + pos0(Re - 1000)*0.44;
}*/


Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::SchillerNaumann::CdRe() const
{
    //volScalarField Re = max(mag(fluid_.U1() - fluid_.U2())*d_/fluid_.thermo1().nu(), 1e-6);

    //return neg(Re_ - 1000)*24.0*(1.0 + 0.15*pow(Re_, 0.687)) + pos0(Re_ - 1000)*0.44*max(Re_, residualRe_);
    return neg(Re_ - 1000)*24.0*(1.0 + 0.15*pow(Re_, 0.687)) + pos0(Re_ - 1000)*0.44*Re_;
}


void Foam::TwoFluidFoam::SchillerNaumann::correct()
{
    //Cd_ = Cd();
    //Ki_ = ((9.0/2.0)*fluid_.thermo1().mu()/(sqr(d_/2)*Cc(d_)))*0.44;//*pos(fluid_.alpha() - 1e-26);
    nu_ = fluid_.thermo1().nu();
    Re_ = mag(fluid_.U1() - fluid_.U2())*d_/nu_;
    Ki_ = 0.75*CdRe()*fluid_.thermo1().rho()*nu_/sqr(d_);
}


Foam::scalar Foam::TwoFluidFoam::SchillerNaumann::dKidp(const label celli) const
{
    /*const volScalarField rho1 = fluid_.thermo1().rho();

    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];

    const scalar drho1dp = rho1[celli]*fluid_.gasProps1().beta_T(fluid_.p()[celli], fluid_.T1()[celli]);

    return 0.75*(Cd_[celli]/d_[celli])*drho1dp*mag(U1 - U2);*/
    return 0.0;
}


Foam::scalar Foam::TwoFluidFoam::SchillerNaumann::dKidalpha(const label celli) const
{
    return 0.0;
}


Foam::vector Foam::TwoFluidFoam::SchillerNaumann::dKidU1(const label celli) const
{
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];
    const scalar dU12 = mag(U1 - U2);
    const scalar d = d_[celli];
    const scalar Re = Re_[celli];
    const scalar nu = nu_[celli];

    const volScalarField& rho1 = fluid_.thermo1().rho();

    //const vector dRedU1 = (d*(U1 - U2))/(nu*mag(U1 - U2));
    //const vector dCdReDU1 = dCdRedx(Re_[celli], dRedU1);

    return 0.75*rho1[celli]*(nu/sqr(d))*dCdRedx(Re, (d*(U1 - U2))/(nu*(dU12 + VSMALL)));

    /*const volScalarField rho1 = fluid_.thermo1().rho();
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];
    return 0.75*(Cd_[celli]/d_[celli])*rho1[celli]*(U1 - U2)/(mag(U1 - U2) + VSMALL);*/
    
    //return vector::zero;
}


Foam::vector Foam::TwoFluidFoam::SchillerNaumann::dKidU2(const label celli) const
{
    return -dKidU1(celli);
}


Foam::scalar Foam::TwoFluidFoam::SchillerNaumann::dKidT1(const label celli) const
{
    /*const volScalarField rho1 = fluid_.thermo1().rho();

    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];

    const scalar drho1dT = -rho1[celli]*fluid_.gasProps1().beta_p(fluid_.p()[celli], fluid_.T1()[celli]);

    return 0.75*(Cd_[celli]/d_[celli])*drho1dT*mag(U1 - U2);*/
    return 0.0;
}


Foam::scalar Foam::TwoFluidFoam::SchillerNaumann::dKidT2(const label celli) const
{
    return 0.0;
}


// ************************************************************************* //
