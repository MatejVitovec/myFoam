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

 // * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
 
bool Foam::functionObjects::totalPressure::calc()
{
    const word thermoName =
    IOobject::groupName
    (
        "thermophysicalProperties",
        phaseName_
    );

    if
    (
        foundObject<volVectorField>(fieldName_)
     && foundObject<fluidThermo>(thermoName)
    )
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(thermoName);

        autoPtr<gasProperties> pGasProps(gasProperties::New(thermo));
        gasProperties& gasProps = pGasProps.ref();
        
        const volVectorField& U = lookupObject<volVectorField>(fieldName_);
        const volScalarField& p = thermo.p();
        const volScalarField& T = thermo.T();

        auto tp0 = tmp<volScalarField>::New
        (
            IOobject
            (
                resultName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", p.dimensions(), 0.0)
        );

        volScalarField& p0 = tp0.ref();

        forAll(p0, celli)
        {
            p0[celli] = this->pTot(p[celli], T[celli], U[celli], gasProps);

            /*const scalar h = gasProps.Hs(p[celli], T[celli]);
            const scalar s = gasProps.S (p[celli], T[celli]);

            const scalar h0 = h + 0.5*magSqr(U[celli]);

            Info << h0 << " , " << h << " , " << s << " , " << p[celli] << endl;

            p0[celli] = gasProps.pHS(h0, s, p[celli]);*/
        }

        forAll(p0.boundaryField(), patchi)
        {
            auto& pp0  = p0.boundaryFieldRef()[patchi];
            const auto& pU = U.boundaryField()[patchi];
            const auto& pp = p.boundaryField()[patchi];
            const auto& pT = T.boundaryField()[patchi];

            forAll(pp0, facei)
            {
                pp0[facei] = this->pTot(pp[facei], pT[facei], pU[facei], gasProps);
            }
        }

        //p0.correctBoundaryConditions();

        return store(resultName_, tp0);
    }

    return false;
}


// * * * * * * * * * * * * Constructor * * * * * * * * * * * * //

Foam::functionObjects::totalPressure::totalPressure
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    phaseName_("")
{
    read(dict);
}


// * * * * * * * * * * * * Read * * * * * * * * * * * * //

bool Foam::functionObjects::totalPressure::read(const dictionary& dict)
{
    fieldExpression::read(dict);

    phaseName_  = dict.lookupOrDefault<word>("phaseName", "");

    resultName_ = IOobject::groupName("p0", phaseName_);

    return true;
}


Foam::scalar Foam::functionObjects::totalPressure::pTot
(
    Foam::scalar p, Foam::scalar T, Foam::vector U, Foam::gasProperties& gasProps
)
{
    scalar S = gasProps.S(p, T);
    scalar H = gasProps.Hs(p, T) + 0.5*magSqr(U);
    scalar T0 = T;
    scalar p0 = p;

    const scalar tol = 1.e-8;
    const label maxIter = 100;
    
    label iter = 0;
    scalar dp, dT;
    do
    {
        scalar Cp = gasProps.Cp(p0, T0);
        scalar beta_p = gasProps.beta_p(p0, T0);
        scalar v = 1.0/gasProps.rho(p0, T0);
        
        scalar dH = H - gasProps.Hs(p0, T0);
        scalar dS = S - gasProps.S(p0, T0);
        dT = dH/Cp;
        dp = (Cp/T0*dT - dS)/(v*beta_p);
        
        
        T0 += dT;
        p0 += dp;
        
        if (iter++ > maxIter)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter
                    << " T  : " << T0
                    << " p  : " << p0
                    << " Z  : " << gasProps.Z(p0,T0)
                    << " Cp : " << gasProps.Cp(p0,T0)
                    << " tol: " << tol
                    << abort(FatalError);
        }
        
    } while ( (mag(dp) > p*tol) || (mag(dT) > T*tol) );
    
    return p0;
}



// ************************************************************************* //
