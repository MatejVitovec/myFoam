#include "totalPressure.H"

#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvc.H"
#include "OFstream.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(totalPressure, 0);
    addToRunTimeSelectionTable(functionObject, totalPressure, dictionary);
}
}


// * * * * * * * * * * * * Constructor * * * * * * * * * * * * //

Foam::functionObjects::totalPressure::totalPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phaseName_("phase"),
    pName_("p"),
    TName_("T"),
    UName_("U"),
    rhoName_("rho"),
    patches_(),
    writeField_(true)
{
    read(dict);
}


// * * * * * * * * * * * * Read * * * * * * * * * * * * //

bool Foam::functionObjects::totalPressure::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    phaseName_  = dict.lookupOrDefault<word>("phase", "phase");
    pName_      = dict.lookupOrDefault<word>("p", "p");
    TName_      = dict.lookupOrDefault<word>("T", "T");
    UName_      = dict.lookupOrDefault<word>("U", "U");

    patches_    = dict.lookupOrDefault<wordList>("patches", wordList());

    writeField_ = dict.lookupOrDefault<bool>("writeField", true);

    return true;
}

// * * * * * * * * * * * Execute * * * * * * * * * * * * //

bool Foam::functionObjects::totalPressure::execute()
{
    const fvMesh& mesh =
        refCast<const fvMesh>(obr_);

    const volScalarField& p =
        mesh.lookupObject<volScalarField>(pName_);

    const volScalarField& T =
        mesh.lookupObject<volScalarField>(TName_);

    const volVectorField& U =
        mesh.lookupObject<volVectorField>(UName_);

    const fluidThermo& thermo =
        lookupObject<fluidThermo>(IOobject::groupName("thermophysicalProperties", phaseName_));

    const autoPtr<gasProperties> pGasProps = gasProperties::New(thermo);
    const gasProperties& gasProps = pGasProps();

    // Create total-pressure field
    volScalarField p0
    (
        IOobject
        (
            "p0",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero",
            p.dimensions(),
            0.0
        )
    );

    // Calculate local total pressure
    forAll(p0, celli)
    {
        const scalar h = gasProps.Hs(p[celli], T[celli]);
        const scalar s = gasProps.S(p[celli], T[celli]);

        const scalar h0 = h + 0.5*magSqr(U[celli]);
        p0[celli] = gasProps.pHS(h0, s, p[celli]);
    }

    p0.correctBoundaryConditions();

    if (writeField_)
    {
        p0.write();
    }

    // Patch statistics
    //calculatePatchAverage(p0);

    return true;
}


// * * * * * * * * * * Patch average * * * * * * * * * * * //

/*void totalPressure::calculatePatchAverage
(
    const volScalarField& p0
) const
{
    const fvMesh& mesh =
        refCast<const fvMesh>(obr_);

    const volVectorField& U =
        mesh.lookupObject<volVectorField>(UName_);

    const volScalarField& rho =
        mesh.lookupObject<volScalarField>(rhoName_);


    forAll(patches_, patchi)
    {
        const word& patchName =
            patches_[patchi];

        const label patchID =
            mesh.boundaryMesh().findPatchID(patchName);


        if (patchID < 0)
        {
            WarningInFunction
                << "Patch " << patchName
                << " was not found."
                << nl;

            continue;
        }


        const vectorField& Sf =
            mesh.Sf().boundaryField()[patchID];

        const fvPatchVectorField& Up =
            U.boundaryField()[patchID];

        const fvPatchScalarField& rhop =
            rho.boundaryField()[patchID];

        const fvPatchScalarField& p0p =
            p0.boundaryField()[patchID];


        scalar massFlow = 0.0;
        scalar p0MassFlow = 0.0;


        forAll(Sf, facei)
        {
            const scalar mdot =
                rhop[facei]
              * (Up[facei] & Sf[facei]);


            massFlow += mdot;

            p0MassFlow +=
                mdot*p0p[facei];
        }


        reduce(massFlow, sumOp<scalar>());
        reduce(p0MassFlow, sumOp<scalar>());


        if (mag(massFlow) > SMALL)
        {
            const scalar p0Mean =
                p0MassFlow/massFlow;


            Info<< "totalPressure:"
                << " patch = " << patchName
                << "  massFlow = " << massFlow
                << "  p0 = " << p0Mean
                << " Pa"
                << nl;


            // Create output directory

            fileName outputDir =
                mesh.time().path()
              / "postProcessing"
              / name();

            mkDir(outputDir);


            fileName outputFile =
                outputDir/(patchName + ".dat");


            OFstream os
            (
                outputFile,
                IOstreamOption::APPEND
            );


            os
                << mesh.time().value()
                << " "
                << p0Mean
                << " "
                << massFlow
                << nl;
        }
    }
}*/


// * * * * * * * * * * * Write * * * * * * * * * * * * //

bool Foam::functionObjects::totalPressure::write()
{
    return true;
}


// ************************************************************************* //
