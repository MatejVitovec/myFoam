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

#include "condensationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "applyFunctions.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationModel, 0);
defineRunTimeSelectionTable(condensationModel, params);


autoPtr<condensationModel>
condensationModel::New
(
    volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const liquidProperties& liquidProps
)
{
    word condensationModelName
    (
        gasThermo.lookup("condensationModel")
    );

    Info<< "Selecting condensation model "
        << condensationModelName << endl;

    auto cstrIter = paramsConstructorTablePtr_->cfind(condensationModelName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown condesnationModel "
            << condensationModelName << endl << endl
            << "Valid condensationModels are : " << endl
            << paramsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<condensationModel>(cstrIter()(alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, liquidProps));
}


condensationModel::condensationModel
(
    volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const liquidProperties& liquidProps
)
:
    mesh_(U.mesh()),
    time_(U.time()),
    alpha_(alpha),
    rho_(rho),
    U_(U),
    alphaRhoPhi_(alphaRhoPhi),
    gasThermo_(gasThermo),
    liquidThermo_(liquidThermo),
    pGasProps_(gasProperties::New(gasThermo)),
    gasProps_(pGasProps_()),
    liquidProps_(liquidProps),
    pSaturation_(saturationCurve::New(gasThermo)),
    saturation_(pSaturation_()),
    nucleationRateMassSource_(
        IOobject
        (
            "mDotNucleation",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimDensity/dimTime, 0.0)
    ),
    growthRateMassSource_(
        IOobject
        (
            "mDotGrowthRate",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimDensity/dimTime, 0.0)
    )
{}


tmp<volScalarField> condensationModel::L() const
{
    volScalarField Ts = saturation_.Ts(liquidThermo_.p());

    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        { return gasProps_.rho(pp,TT); }; //TODO namam gasProps
        
    volScalarField rhos = applyFunction2(f, liquidThermo_.p(), Ts, "rhos", dimDensity);

    return tmp<volScalarField>
        (
	        //new volScalarField("L", saturation_.dpsdT(Ts)*Ts/rhos)
            new volScalarField("L", gasThermo_.p()*(gasThermo_.T() - Ts) + saturation_.dpsdT(Ts)*Ts/rhos)
        );

    /*const std::function<scalar(scalar)> f = [&](scalar TT)
        {return 461.52*(-2.7246e-2*sqr(TT) + 2*1.6853e-5*pow(TT,3) + 2.4576*TT + 6094.4642);};

    return applyFunction1(f, T(), "L", dimEnergy/dimMass);*/
}


tmp<volScalarField> condensationModel::w() const
{
    return (liquidThermo_.rho()*(1.0 - alpha_))/(gasThermo_.rho() + (1.0 - alpha_)*(liquidThermo_.rho() - gasThermo_.rho()));
}

tmp<volScalarField> condensationModel::dropletDiameter() const
{
    return
        tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "condensationDropletDiemeter",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                2.0
            )
        );
}

}
}
// ************************************************************************* //
