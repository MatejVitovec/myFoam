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

#include "nucleationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{
    defineTypeNameAndDebug(nucleationModel, 0);
    defineRunTimeSelectionTable(nucleationModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WetSteam::nucleationModel::nucleationModel
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur,
    const gasProperties& gasProps,
    const liquidProperties& liquidProps
)
:
    gasThermo_(gasThermo),
    liquidThermo_(liquidThermo),
    saturation_(satur),
    gasProps_(gasProps),
    liquidProps_(liquidProps),
    rc_(
        IOobject
        (
            "rc",
            gasThermo_.p().mesh().time().timeName(),
            gasThermo_.p().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        gasThermo_.p().mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    J_(
        IOobject
        (
            "J",
            gasThermo_.p().mesh().time().timeName(),
            gasThermo_.p().mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        gasThermo_.p().mesh(),
        dimensionedScalar("zero", dimless/dimTime/dimVolume, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::WetSteam::nucleationModel>
Foam::WetSteam::nucleationModel::New
(
    const dictionary& dict,
    const fluidThermo& gasThermo,
    const fluidThermo& liquidThermo,
    const saturation& satur,
    const gasProperties& gasProps,
    const liquidProperties& liquidProps
)
{
    word nucleationModelName(dict.lookup("nucleationModel"));

    Info<< "Selecting nucleationModel "
        << nucleationModelName << endl;

    auto* ctorPtr = dictionaryConstructorTable(nucleationModelName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "nucleationModel",
            nucleationModelName,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(dict, gasThermo, liquidThermo, satur, gasProps, liquidProps);
}

// * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::WetSteam::nucleationModel::sigma() const
{
    const fvMesh& mesh = gasThermo_.p().mesh();
    const volScalarField& T_g = gasThermo_.T();

    tmp<volScalarField> tsigma
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("sigmaZero", dimMass/dimTime/dimTime, 0.0)
        )
    );

    volScalarField& sigma = tsigma.ref();

    forAll(sigma, i)
    {
        const scalar tau = max(1.0 - T_g[i]/647.096, 0.0);
        sigma[i] = 235.8e-3*pow(tau, 1.256)*(1 - 0.625*tau);
        //sigma[i] = liquidProps_.sigma(p[i], T_g[i]);
    }

    return tsigma;
}

Foam::tmp<Foam::volScalarField>
Foam::WetSteam::nucleationModel::criticalDropletRadius(const volScalarField& sigma) const
{
    const fvMesh& mesh = gasThermo_.p().mesh();

    tmp<volScalarField> trc
    (
        new volScalarField 
        (
            IOobject
            (
                "rc",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimLength, 0.0)
        )
    );

    volScalarField& rc = trc.ref();

    const volScalarField& p = gasThermo_.p();
    const volScalarField& T_g = gasThermo_.T();
    const volScalarField& ps = saturation_.ps();
    //const volScalarField& rho_l = liquidThermo_.rho();
    const volScalarField rhos_l = saturation_.rhosl();

    forAll(rc, i)
    {
        const scalar dH = gasProps_.Ha(p[i], T_g[i]) - gasProps_.Ha(ps[i], T_g[i]);
        const scalar ds = gasProps_.S( p[i], T_g[i]) - gasProps_.S( ps[i], T_g[i]);
        const scalar dG = dH - T_g[i]*ds;
        rc[i] = max(2*sigma[i]/(rhos_l[i]*dG), 0.0);
        //rc[i] = 2*sigma[i]/(rho_l[i]*dG);

        //rc[i] = 2*sigma/(L*rho_l[i]*log(Ts[i]/T_g[i]));

        //const scalar S = p[i]/ps[i];
        //rc[i] = 2*sigma/(rho_l[i]*Rg*T_g[i]*log(S));
    }

    return trc;
}

const Foam::volScalarField&
Foam::WetSteam::nucleationModel::J() const
{
    return J_;
}

const Foam::volScalarField&
Foam::WetSteam::nucleationModel::rc() const
{
    return rc_;
}

// ************************************************************************* //
