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
    dbFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "twoFluid.H"
#include "twoFluidConvectiveFlux.H"
#include "twoFluidFvc.H"
#include "dragModel.H"

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

        scalar initialRezAlphaRho1 = 0, initialRezAlphaRho2 = 0, initialRezAlphaRhoU1 = 0, initialRezAlphaRhoU2 = 0, initialRezEpsilon1 = 0, initialRezEpsilon2 = 0;


        for (int intIter = 0; intIter < lusgsIntIters; intIter++)
        {
            twoFluidFlux.computeFlux();

            volScalarField rho = alpha1*rho1 + alpha2*rho2;
            volVectorField virtualVelocity = (alpha1*rho1*U1 + alpha2*rho2*U2)/rho;
            volVectorField dragTerm = drag.K(d)*(fluid.U1() - fluid.U2());

            volVectorField virtualMassTerm= -0.5*rho*alpha1*alpha2*(DU2 - DU1);

            volScalarField dAlphaRho1(-dt*(
                fvc::ddt(conservative.alphaRho1()) + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg())
            ));

            volScalarField dAlphaRho2(-dt*(
                fvc::ddt(conservative.alphaRho2()) + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg())
            ));

            volVectorField dAlphaRhoU1(-dt*(
                fvc::ddt(conservative.alphaRhoU1())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha1_pos()*mesh.Sf(), twoFluidFlux.alpha1_neg()*mesh.Sf())
                //+ drag.K(d)*(fluid.U1() - fluid.U2())
                + dragTerm
                //+ virtualMassTerm
            ));

            volVectorField dAlphaRhoU2(-dt*(
                fvc::ddt(conservative.alphaRhoU2())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha2_pos()*mesh.Sf(), twoFluidFlux.alpha2_neg()*mesh.Sf())
                //+ drag.K(d)*(fluid.U2() - fluid.U1())
                - dragTerm
                //- virtualMassTerm
            ));

            volScalarField dEpsilon1(-dt*(
                fvc::ddt(conservative.epsilon1())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg())
                //+ (dragTerm & U1)
                + ((dragTerm /*+ virtualMassTerm*/) & virtualVelocity)
            ));

            volScalarField dEpsilon2(-dt*(
                fvc::ddt(conservative.epsilon2())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg())
                //- (dragTerm & U2)
                - ((dragTerm /*+ virtualMassTerm*/) & virtualVelocity)
            ));

            scalar rezAlphaRho1  = fvc::domainIntegrate(mag(dAlphaRho1) /dt).value();
            scalar rezAlphaRho2  = fvc::domainIntegrate(mag(dAlphaRho2) /dt).value();
            scalar rezAlphaRhoU1 = fvc::domainIntegrate(mag(dAlphaRhoU1)/dt).value();
            scalar rezAlphaRhoU2 = fvc::domainIntegrate(mag(dAlphaRhoU2)/dt).value();
            scalar rezEpsilon1   = fvc::domainIntegrate(mag(dEpsilon1)  /dt).value();
            scalar rezEpsilon2   = fvc::domainIntegrate(mag(dEpsilon2)  /dt).value();

            if (intIter == 0)
            {
                initialRezAlphaRho1  = rezAlphaRho1;
                initialRezAlphaRho2  = rezAlphaRho2;
                initialRezAlphaRhoU1 = rezAlphaRhoU1;
                initialRezAlphaRhoU2 = rezAlphaRhoU2;
                initialRezEpsilon1   = rezEpsilon1;
                initialRezEpsilon2   = rezEpsilon2;
            }

#           include "lusgsSweep.H"

            conservative.alphaRho1()  += dAlphaRho1;
            conservative.alphaRho2()  += dAlphaRho2;
            conservative.alphaRhoU1() += dAlphaRhoU1;
            conservative.alphaRhoU2() += dAlphaRhoU2;
            conservative.epsilon1()   += dEpsilon1;
            conservative.epsilon2()   += dEpsilon2;


            fluid.correct();
            fluid.blendVanishingFluid();
            //fluid.bound();
            fluid.correctBoundaryCondition();

            /*if (runTime.timeIndex() > 12570)
            {
                Info << ">>> Forcing write , alpha: <<<" << alpha[15248] << endl;
                runTime.writeNow();
                runTime.write();
            }*/

            fluid.correctThermo();
            fluid.correctInterfacialPressure();
            fluid.correctConservative();

            DU1 = fvc::ddt(U1) + fvc::div(fvc::flux(U1), U1);
            DU2 = fvc::ddt(U2) + fvc::div(fvc::flux(U2), U2);

            scalar finalRezAlphaRho1  = fvc::domainIntegrate(mag(dAlphaRho1) /dt).value();
            scalar finalRezAlphaRho2  = fvc::domainIntegrate(mag(dAlphaRho2) /dt).value();
            scalar finalRezAlphaRhoU1 = fvc::domainIntegrate(mag(dAlphaRhoU1)/dt).value();
            scalar finalRezAlphaRhoU2 = fvc::domainIntegrate(mag(dAlphaRhoU2)/dt).value();
            scalar finalRezEpsilon1   = fvc::domainIntegrate(mag(dEpsilon1)  /dt).value();
            scalar finalRezEpsilon2   = fvc::domainIntegrate(mag(dEpsilon2)  /dt).value();

            Info << "LUSGS:  Solving for alphaRho1,          " << "Initial residual = " << rezAlphaRho1  << ", " << "Final residual = " << finalRezAlphaRho1  << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for alphaRho2,          " << "Initial residual = " << rezAlphaRho2  << ", " << "Final residual = " << finalRezAlphaRho2  << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for alphaRhoU1,         " << "Initial residual = " << rezAlphaRhoU1 << ", " << "Final residual = " << finalRezAlphaRhoU1 << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for alphaRhoU2,         " << "Initial residual = " << rezAlphaRhoU2 << ", " << "Final residual = " << finalRezAlphaRhoU2 << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for alpha1(rhoE1+pInt), " << "Initial residual = " << rezEpsilon1   << ", " << "Final residual = " << finalRezEpsilon1   << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for alpha2(rhoE2+pInt), " << "Initial residual = " << rezEpsilon2   << ", " << "Final residual = " << finalRezEpsilon2   << ", No Iterations 1" << nl;

            bool lastIteration = ((intIter + 1) == lusgsIntIters);

            lastIteration = lastIteration || ( 
                (finalRezAlphaRho1  < lusgsRelTol * initialRezAlphaRho1) && 
                (finalRezAlphaRho2  < lusgsRelTol * initialRezAlphaRho2) && 
                (finalRezAlphaRhoU1 < lusgsRelTol * initialRezAlphaRhoU1) && 
                (finalRezAlphaRhoU2 < lusgsRelTol * initialRezAlphaRhoU2) && 
                (finalRezEpsilon1   < lusgsRelTol * initialRezEpsilon1) && 
                (finalRezEpsilon2   < lusgsRelTol * initialRezEpsilon2) );

            lastIteration = lastIteration || ( 
                (finalRezAlphaRho1  < lusgsTolerance) && 
                (finalRezAlphaRho2  < lusgsTolerance) && 
                (finalRezAlphaRhoU1 < lusgsTolerance) && 
                (finalRezAlphaRhoU2 < lusgsTolerance) && 
                (finalRezEpsilon1   < lusgsTolerance) && 
                (finalRezEpsilon2   < lusgsTolerance) );

            if (lastIteration) {
                Info << "LUSGS: converged in " << intIter + 1 << " iterations." << nl << nl;
                break;
            }
        }

        rho1.ref() = thermo1.rho();
        rho2.ref() = thermo2.rho();

        runTime.write();

        //runTime.printExecutionTime(Info);

        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
