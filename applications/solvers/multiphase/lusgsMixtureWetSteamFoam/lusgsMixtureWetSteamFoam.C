/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

Application
    lusgsTwoFluidFoam

Group
    grpCompressibleSolvers

Description
    Density-based lusgs mixture wet-steam flow solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvcSmooth.H"

#include "compressibleMixture.H"
#include "compressibleMixtureConvectiveFlux.H"
#include "condensationModel.H"
#include "saturation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " Riemann problem flux calculation"
    );

    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readLUSGSControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        amaxSf = fluid.amaxSf();
        #include "courantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << "Iter = " << runTime.timeIndex() << endl;
        Info << "Time = " << runTime.timeName() << endl;
        
        dimensionedScalar dt = runTime.deltaT();

        //for (int intIter = 0; intIter < lusgsIntIters; intIter++)
        {
            mixtureFlux.computeFlux();

            condensationMassSource = condensation.condensationRateMassSource();

            volScalarField resRho (-dt*(fvc::div(mixtureFlux.rhoFlux())));
            volVectorField resRhoU(-dt*(fvc::div(mixtureFlux.rhoUFlux())));
            volScalarField resRhoE(-dt*(fvc::div(mixtureFlux.rhoEFlux())));
            volScalarField resRhow(-dt*(fvc::div(mixtureFlux.rhowFlux()) - condensationMassSource));

            dp  = dimensionedScalar("pzero", dimPressure, 0.0);
            dU  = dimensionedVector("Uzero", dimVelocity, vector(0, 0, 0));
            dT  = dimensionedScalar("Tzero", dimTemperature, 0.0);
            dw  = dimensionedScalar("wzero", dimless, 0.0);

            #include "lusgsSweep.H"

            p += dp;
            U += dU;
            T += dT;
            w += dw;

            //Info << dp << endl;
            //std::cin.ignore();

            //fluid.correct();
            //fluid.bound();
            fluid.constrainW();
            fluid.correctBoundaryCondition();
            fluid.correctThermo();
            //fluid.correctConservative();

            rho   = fluid.rho();
            alpha = fluid.alpha();

            //Info << "alpha min: " << gMin(alpha) << " alpha max: " << gMax(alpha) << " w min: " <<gMin(fluid.w()) << " w max: " << gMax(fluid.w()) << endl;

            satur.correct();

            scalar finalRezp    = fvc::domainIntegrate(mag(dp)/dt).value();
            scalar finalRezU    = fvc::domainIntegrate(mag(dU)/dt).value();
            scalar finalRezT    = fvc::domainIntegrate(mag(dT)/dt).value();
            scalar finalRezw    = fvc::domainIntegrate(mag(dw)/dt).value();

            Info << "LUSGS:  Solving for p, " << "Final residual = " << finalRezp   << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for U, " << "Final residual = " << finalRezU   << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for T, " << "Final residual = " << finalRezT   << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for w, " << "Final residual = " << finalRezw   << ", No Iterations " << lusgsIntIters << nl;
        }

        condensation.correct();

        runTime.write();

        if (mesh.time().outputTime())
        {
            rhov.write();
            rhol.write();
            rho.write();
            //e1.write();
            //e2.write();
        }

        /*if (runTime.timeIndex() > 21650)
        {
            Info << ">>> Forcing write <<<" << endl;
            runTime.writeNow();
            runTime.write();
        }*/

        //runTime.printExecutionTime(Info);

        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
