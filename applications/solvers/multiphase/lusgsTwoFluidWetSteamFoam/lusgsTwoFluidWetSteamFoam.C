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
    Density-based lusgs two-fluid compressible flow solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "twoFluid.H"
#include "twoFluidConvectiveFlux.H"
#include "twoFluidFvc.H"
#include "dragModel.H"
#include "condensationModel.H"

#include <eigen3/Eigen/Dense>

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

        volVectorField dragSource
        (
            "dragS",
            drag.K(condensation.dropletDiameter())*(U1 - U2)
        );

        //for (int intIter = 0; intIter < lusgsIntIters; intIter++)
        {
            twoFluidFlux.computeFlux();

            volVectorField Uint = (alpha1*rho1*U1 + alpha2*rho2*U2)/(alpha1*rho1 + alpha2*rho2);

            volScalarField Ts = condensation.saturation().Ts(p);
            volScalarField Hvint = condensation.saturation().hsv(Ts) + (Uint & U1) - 0.5*magSqr(U1);
            volScalarField Hlint = condensation.saturation().hsl(Ts) + (Uint & U2) - 0.5*magSqr(U2);

            volScalarField dragK = drag.K(condensation.dropletDiameter());
            dragSource = drag.K(condensation.dropletDiameter())*(U1 - U2);

            volScalarField condensationMassSource = condensation.condensationRateMassSource();
            volVectorField condensationMomentumSource = Uint*condensation.condensationRateMassSource();
            volScalarField condensationEnergyVaporSource = condensation.condensationRateMassSource()*(Hvint - condensation.L());
            volScalarField condensationEnergyLiquidSource = condensation.condensationRateMassSource()*Hlint;

            volScalarField rezAlphaRho1(-dt*(
                /*fvc::ddt(conservative.alphaRho1())
                +*/ TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg())
                + condensationMassSource
            ));

            volScalarField rezAlphaRho2(-dt*(
                /*fvc::ddt(conservative.alphaRho2())
                +*/ TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg())
                - condensationMassSource
            ));

            volVectorField rezAlphaRhoU1(-dt*(
                /*fvc::ddt(conservative.alphaRhoU1())
                +*/ TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha1_pos()*mesh.Sf(), twoFluidFlux.alpha1_neg()*mesh.Sf())
                + dragSource
                + condensationMomentumSource
            ));

            volVectorField rezAlphaRhoU2(-dt*(
                /*fvc::ddt(conservative.alphaRhoU2())
                +*/ TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha2_pos()*mesh.Sf(), twoFluidFlux.alpha2_neg()*mesh.Sf())
                - dragSource
                - condensationMomentumSource
            ));

            volScalarField rezEpsilon1(-dt*(
                /*fvc::ddt(conservative.epsilon1())
                +*/ TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg())
                //+ (dragSource & Uint)
                + (dragSource & U1)
                + condensationEnergyVaporSource                
            ));

            volScalarField rezEpsilon2(-dt*(
                /*fvc::ddt(conservative.epsilon2())
                +*/ TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg())
                //- (dragSource & Uint)
                - (dragSource & U2)
                - condensationEnergyLiquidSource
            ));

            dp     = dimensionedScalar("pzero", dimPressure, 0.0);
            dalpha = dimensionedScalar("alphazero", dimless, 0.0);
            dU1    = dimensionedVector("U1zero", dimVelocity, vector(0, 0, 0));
            dU2    = dimensionedVector("U2zero", dimVelocity, vector(0, 0, 0));
            dT1    = dimensionedScalar("T1zero", dimTemperature, 0.0);
            dT2    = dimensionedScalar("T2zero", dimTemperature, 0.0);

            #include "lusgsSweep.H"

            p     += dp;
            alpha += dalpha;
            U1    += dU1;
            U2    += dU2;
            T1    += dT1;
            T2    += dT2;

            //fluid.correct();
            //fluid.blendVanishingFluid();
            fluid.blendVanishingFluid(Ts);
            //fluid.bound();
            fluid.correctBoundaryCondition();
            fluid.correctThermo();
            fluid.correctInterfacialPressure();
            //fluid.correctConservative();

            scalar finalRezp     = fvc::domainIntegrate(mag(dp)    /dt).value();
            scalar finalRezalpha = fvc::domainIntegrate(mag(dalpha)/dt).value();
            scalar finalRezU1    = fvc::domainIntegrate(mag(dU1)   /dt).value();
            scalar finalRezU2    = fvc::domainIntegrate(mag(dU2)   /dt).value();
            scalar finalRezT1    = fvc::domainIntegrate(mag(dT1)   /dt).value();
            scalar finalRezT2    = fvc::domainIntegrate(mag(dT2)   /dt).value();

            Info << "LUSGS:  Solving for p,     " << "Final residual = " << finalRezp     << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for alpha, " << "Final residual = " << finalRezalpha << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for U1,    " << "Final residual = " << finalRezU1    << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for U2,    " << "Final residual = " << finalRezU2    << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for T1,    " << "Final residual = " << finalRezT1    << ", No Iterations " << lusgsIntIters << nl;
            Info << "LUSGS:  Solving for T2,    " << "Final residual = " << finalRezT2    << ", No Iterations " << lusgsIntIters << nl;
        }

        rho = fluid.rho();
        rhoMPhi2 = linearInterpolate(rho*U2) & mesh.Sf();
        condensation.correct();


        forAll(s1, celli)
        {
            s1[celli] = fluid.gasProps1().S(p[celli], T1[celli]);
        }
        h1 = e1 + p/rho1;

        runTime.write();

        if (mesh.time().outputTime())
        {
            rho1.write();
            rho2.write();
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
