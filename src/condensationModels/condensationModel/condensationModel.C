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
    const saturation& satur,
    const dictionary& dict
)
{
    const dictionary& condensationDict = dict.subDict("condensation");

    word condensationModelName
    (
        condensationDict.lookup("condensationModel")
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

    return autoPtr<condensationModel>(cstrIter()(alpha, rho, U, alphaRhoPhi, gasThermo, liquidThermo, satur, condensationDict));
}


condensationModel::condensationModel
(
    volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur,
    const dictionary& dict
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
    pliquidProps_(liquidProperties::New(dict.lookupOrDefault<word>("liquidProperties", "H2O"))),
    liquidProps_(pliquidProps_()),
    //pSaturation_(saturationCurve::New(dict)),
    //saturation_(pSaturation_()),
    saturation_(satur),
    Kn_(
        IOobject
        (
            "Kn",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimless, 1.0)
    ),
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


//tmp<volScalarField> condensationModel::L() const
//{
    /*const volScalarField Ts = saturation_.Ts(gasThermo_.p());
    const volScalarField p = gasThermo_.p();

    labelList cells = identity(mesh_.nCells());
    scalarField pCells = gasThermo_.p().internalField();
    scalarField TsCells = Ts.internalField();
    scalarField rhosgCells = gasThermo_.rhoEoS(pCells, TsCells, cells);
    scalarField rhoslCells = liquidThermo_.rhoEoS(pCells, TsCells, cells);
    
    volScalarField rhosg
    (
        IOobject
        (
            "rhosg",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        gasThermo_.rho()
    );

    volScalarField rhosl
    (
        IOobject
        (
            "rhosl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        liquidThermo_.rho()
    );

    forAll(Ts, i)
    {
        rhosg[i] = rhosgCells[i];
        rhosl[i] = rhoslCells[i];
    }*/
   
    /*TADY const volScalarField p = gasThermo_.p();

    return tmp<volScalarField>
        (
            //new volScalarField("L", gasThermo_.Cp()*(gasThermo_.T() - Ts) + (1/rhosg - 1/rhosl)*Ts*saturation_.dpsdT(Ts))
            new volScalarField("L", (gasThermo_.he() + p/gasThermo_.rho()) - saturation_.hsl())
        );*/


    /*const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        { return gasProps_.rho(pp,TT); };
        
    volScalarField rhos = applyFunction2(f, liquidThermo_.p(), Ts, "rhos", dimDensity);

    return tmp<volScalarField>
        (
            //new volScalarField("L", saturation_.dpsdT(Ts)*Ts/rhos)
            new volScalarField("L", gasThermo_.Cp()*(gasThermo_.T() - Ts) + saturation_.dpsdT(Ts)*Ts/rhos)
        );*/

    /////END

    /*const std::function<scalar(scalar)> f = [&](scalar TT)
        {return 461.52*(-2.7246e-2*sqr(TT) + 2*1.6853e-5*pow(TT,3) + 2.4576*TT + 6094.4642);};

    return applyFunction1(f, T(), "L", dimEnergy/dimMass);*/
//}


tmp<volScalarField> condensationModel::w() const
{
    return (liquidThermo_.rho()*alpha_)/(gasThermo_.rho() + alpha_*(liquidThermo_.rho() - gasThermo_.rho()));
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
