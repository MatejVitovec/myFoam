/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "dragModel.H"
#include "twoFluid.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace TwoFluidFoam
{
    defineTypeNameAndDebug(dragModel, 0);
    defineRunTimeSelectionTable(dragModel, dictionary);
}
}

const Foam::dimensionSet Foam::TwoFluidFoam::dragModel::dimK(1, -3, -1, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModel::dragModel
(
    const twoFluid& fluid,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, twoFluid.name()),
            typeName,
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    fluid_(fluid)
{}


Foam::TwoFluidFoam::dragModel::dragModel
(
    const dictionary& dict,
    const twoFluid& fluid,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            //IOobject::groupName(typeName, fluid.name()),
            typeName,
            fluid.mesh().time().timeName(),
            fluid.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    fluid_(fluid)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
Foam::autoPtr<Foam::TwoFluidFoam::dragModel>
Foam::TwoFluidFoam::dragModel::New
(
    const twoFluid& fluid
)
{
    const dictionary& dict = fluid.subDict("drag");

    const word modelType(dict.get<word>("type"));

    Info<< "Selecting dragModel "
        << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "dragModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, fluid, true);
}


Foam::autoPtr<Foam::TwoFluidFoam::dragModel>
Foam::TwoFluidFoam::dragModel::New
(
    const dictionary& dict,
    const twoFluid& fluid
)
{
    const word modelType(dict.get<word>("type"));

    Info<< "Selecting dragModel for "
        << fluid << ": " << modelType << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "dragModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, fluid, true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::TwoFluidFoam::dragModel::~dragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::dragModel::Ki(const volScalarField& d) const
{
    return (0.75*CdRe()*fluid_.thermo1().rho()*mag(fluid_.U1() - fluid_.U2()))/d;
}


Foam::tmp<Foam::volScalarField> Foam::TwoFluidFoam::dragModel::K(const volScalarField& d) const
{
    return Ki(d)*fluid_.alpha();
}





Foam::scalar Foam::TwoFluidFoam::dragModel::dKdp(const label celli, const scalar d, const scalar Cd) const
{
    const volScalarField rho1 = fluid_.thermo1().rho();

    const scalar alpha = fluid_.alpha()[celli];
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];

    const scalar drho1dp = rho1[celli]*fluid_.gasProps1().beta_T(fluid_.p()[celli], fluid_.T1()[celli]);

    return 0.75*Cd*(alpha/d)*drho1dp*mag(U1 - U2);
}

Foam::scalar Foam::TwoFluidFoam::dragModel::dKdalpha(const label celli, const scalar d, const scalar Cd) const
{
    const volScalarField rho1 = fluid_.thermo1().rho();
    
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];

    return 0.75*Cd*(rho1[celli]/d)*mag(U1 - U2);
}

Foam::vector Foam::TwoFluidFoam::dragModel::dKdU1(const label celli, const scalar d, const scalar Cd) const
{
    const volScalarField rho1 = fluid_.thermo1().rho();

    const scalar alpha = fluid_.alpha()[celli];
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];

    return 0.75*(Cd/d)*alpha*rho1[celli]*(U1 - U2)/(mag(U1 - U2) + VSMALL);
}

Foam::vector Foam::TwoFluidFoam::dragModel::dKdU2(const label celli, const scalar d, const scalar Cd) const
{
    return -dKdU1(celli, d, Cd);
}

Foam::scalar Foam::TwoFluidFoam::dragModel::dKdT1(const label celli, const scalar d, const scalar Cd) const
{
    const volScalarField rho1 = fluid_.thermo1().rho();

    const scalar alpha = fluid_.alpha()[celli];
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];

    const scalar drho1dT = -rho1[celli]*fluid_.gasProps1().beta_p(fluid_.p()[celli], fluid_.T1()[celli]);

    return 0.75*Cd*(alpha/d)*drho1dT*mag(U1 - U2);
}

Foam::scalar Foam::TwoFluidFoam::dragModel::dKdT2(const label celli, const scalar d, const scalar Cd) const
{
    return 0.0;
}


std::array<Foam::scalar, 100> Foam::TwoFluidFoam::dragModel::dSdpUT(const label celli, const scalar K, const scalar d) const
{
    constexpr auto idx = [](std::size_t i, std::size_t j) constexpr
    {
        return j * 10 + i;      // column-major
    };

    std::array<Foam::scalar, 100> out{};

    const scalar alpha = fluid_.alpha()[celli];
    const scalar p = fluid_.p()[celli];
    const vector U1 = fluid_.U1()[celli];
    const vector U2 = fluid_.U2()[celli];
    const scalar T1 = fluid_.T1()[celli];
    const scalar T2 = fluid_.T2()[celli];

    const vector U1mU2 = U1 - U2;
    const scalar U1mU2dotU1 = U1mU2 & U1;
    const scalar U2mU1dotU2 = (-U1mU2) & U2;


    const volScalarField CdRe_ = CdRe();
    const scalar Cd = CdRe_[celli];

    const scalar dKdp_ = dKdp(celli, d, Cd);
    const scalar dKdalpha_ = dKdalpha(celli, d, Cd);
    const vector dKdU1_ = dKdU1(celli, d, Cd);
    const vector dKdU2_ = dKdU2(celli, d, Cd);
    const scalar dKdT1_ = dKdT1(celli, d, Cd);
    const scalar dKdT2_ = dKdT2(celli, d, Cd);

    //const scalar dKdp_ = 0.0;
    //const scalar dKdalpha_ = 0.0;
    //const vector dKdU1_ = vector::zero;
    //const vector dKdU2_ = vector::zero;
    //const scalar dKdT1_ = 0.0;
    //const scalar dKdT2_ = 0.0;


    out[idx(1, 0)] =  dKdp_*U1mU2.x();
    out[idx(2, 0)] =  dKdp_*U1mU2.y();
    out[idx(3, 0)] =  dKdp_*U1mU2.z();
    out[idx(4, 0)] =  dKdp_*U1mU2dotU1;
    out[idx(6, 0)] = -dKdp_*U1mU2.x();
    out[idx(7, 0)] = -dKdp_*U1mU2.y();
    out[idx(8, 0)] = -dKdp_*U1mU2.z();
    out[idx(9, 0)] =  dKdp_*U2mU1dotU2;


    out[idx(1, 1)] =  dKdalpha_*U1mU2.x();
    out[idx(2, 1)] =  dKdalpha_*U1mU2.y();
    out[idx(3, 1)] =  dKdalpha_*U1mU2.z();
    out[idx(4, 1)] =  dKdalpha_*U1mU2dotU1;
    out[idx(6, 1)] = -dKdalpha_*U1mU2.x();
    out[idx(7, 1)] = -dKdalpha_*U1mU2.y();
    out[idx(8, 1)] = -dKdalpha_*U1mU2.z();
    out[idx(9, 1)] =  dKdalpha_*U2mU1dotU2;


    out[idx(1, 2)] =  dKdU1_.x()*U1mU2.x()  + K;
    out[idx(2, 2)] =  dKdU1_.x()*U1mU2.y(); //
    out[idx(3, 2)] =  dKdU1_.x()*U1mU2.z(); //
    out[idx(4, 2)] =  dKdU1_.x()*U1mU2dotU1 + K*(2*U1.x() - U2.x());
    out[idx(6, 2)] = -dKdU1_.x()*U1mU2.x()  - K;
    out[idx(7, 2)] = -dKdU1_.x()*U1mU2.y(); //
    out[idx(8, 2)] = -dKdU1_.x()*U1mU2.z(); //
    out[idx(9, 2)] =  dKdU1_.x()*U2mU1dotU2 - K*U2.x();

    out[idx(1, 3)] =  dKdU1_.y()*U1mU2.x(); //o
    out[idx(2, 3)] =  dKdU1_.y()*U1mU2.y()  + K;
    out[idx(3, 3)] =  dKdU1_.y()*U1mU2.z(); //
    out[idx(4, 3)] =  dKdU1_.y()*U1mU2dotU1 + K*(2*U1.y() - U2.y());
    out[idx(6, 3)] = -dKdU1_.y()*U1mU2.x(); //
    out[idx(7, 3)] = -dKdU1_.y()*U1mU2.y()  - K;
    out[idx(8, 3)] = -dKdU1_.y()*U1mU2.z(); //
    out[idx(9, 3)] =  dKdU1_.y()*U2mU1dotU2 - K*U2.y();

    out[idx(1, 4)] =  dKdU1_.z()*U1mU2.x(); //
    out[idx(2, 4)] =  dKdU1_.z()*U1mU2.y(); //
    out[idx(3, 4)] =  dKdU1_.z()*U1mU2.z()  + K;
    out[idx(4, 4)] =  dKdU1_.z()*U1mU2dotU1 + K*(2*U1.z() - U2.z());
    out[idx(6, 4)] = -dKdU1_.z()*U1mU2.x(); //
    out[idx(7, 4)] = -dKdU1_.z()*U1mU2.y(); //
    out[idx(8, 4)] = -dKdU1_.z()*U1mU2.z()  - K;
    out[idx(9, 4)] =  dKdU1_.z()*U2mU1dotU2 - K*U2.z();


    out[idx(1, 5)] =  dKdU2_.x()*U1mU2.x()  - K;
    out[idx(2, 5)] =  dKdU2_.x()*U1mU2.y(); //
    out[idx(3, 5)] =  dKdU2_.x()*U1mU2.z(); //
    out[idx(4, 5)] =  dKdU2_.x()*U1mU2dotU1 - K*U1.x(); 
    out[idx(6, 5)] = -dKdU2_.x()*U1mU2.x()  + K;
    out[idx(7, 5)] = -dKdU2_.x()*U1mU2.y(); //
    out[idx(8, 5)] = -dKdU2_.x()*U1mU2.z(); //
    out[idx(9, 5)] =  dKdU2_.x()*U2mU1dotU2  + K*(2*U2.x() - U1.x());

    out[idx(1, 6)] =  dKdU2_.y()*U1mU2.x(); //
    out[idx(2, 6)] =  dKdU2_.y()*U1mU2.y()  - K;
    out[idx(3, 6)] =  dKdU2_.y()*U1mU2.z(); //
    out[idx(4, 6)] =  dKdU2_.y()*U1mU2dotU1 - K*U1.y(); 
    out[idx(6, 6)] = -dKdU2_.y()*U1mU2.x(); //
    out[idx(7, 6)] = -dKdU2_.y()*U1mU2.y()  + K;
    out[idx(8, 6)] = -dKdU2_.y()*U1mU2.z(); //
    out[idx(9, 6)] =  dKdU2_.y()*U2mU1dotU2 + K*(2*U2.y() - U1.y());

    out[idx(1, 7)] =  dKdU2_.z()*U1mU2.x(); //
    out[idx(2, 7)] =  dKdU2_.z()*U1mU2.y(); //
    out[idx(3, 7)] =  dKdU2_.z()*U1mU2.z()  - K;
    out[idx(4, 7)] =  dKdU2_.z()*U1mU2dotU1 - K*U1.z(); 
    out[idx(6, 7)] = -dKdU2_.z()*U1mU2.x(); //
    out[idx(7, 7)] = -dKdU2_.z()*U1mU2.y(); //
    out[idx(8, 7)] = -dKdU2_.z()*U1mU2.z()  + K;
    out[idx(9, 7)] =  dKdU2_.z()*U2mU1dotU2 + K*(2*U2.z() - U1.z());


    out[idx(1, 8)] =  dKdT1_*U1mU2.x();
    out[idx(2, 8)] =  dKdT1_*U1mU2.y();
    out[idx(3, 8)] =  dKdT1_*U1mU2.z();
    out[idx(4, 8)] =  dKdT1_*U1mU2dotU1;
    out[idx(6, 8)] = -dKdT1_*U1mU2.x();
    out[idx(7, 8)] = -dKdT1_*U1mU2.y();
    out[idx(8, 8)] = -dKdT1_*U1mU2.z();
    out[idx(9, 8)] =  dKdT1_*U2mU1dotU2;

    out[idx(1, 9)] =  dKdT2_*U1mU2.x();
    out[idx(2, 9)] =  dKdT2_*U1mU2.y();
    out[idx(3, 9)] =  dKdT2_*U1mU2.z();
    out[idx(4, 9)] =  dKdT2_*U1mU2dotU1;
    out[idx(6, 9)] = -dKdT2_*U1mU2.x();
    out[idx(7, 9)] = -dKdT2_*U1mU2.y();
    out[idx(8, 9)] = -dKdT2_*U1mU2.z();
    out[idx(9, 9)] =  dKdT2_*U2mU1dotU2;

    return out;
}


bool Foam::TwoFluidFoam::dragModel::writeData(Ostream& os) const
{
    return os.good();
}


// ************************************************************************* //
