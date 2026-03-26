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

#include <eigen3/Eigen/Dense>
#include "emptyPolyPatch.H"

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
    #include "createLUSGSControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readLUSGSControls.H"

        amaxSf = fluid.amaxSf();
        //amaxSf = mag(fvc::interpolate(U1) & mesh.Sf()) + mesh.magSf()*fvc::interpolate(a1);
        #include "courantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << "Iter = " << runTime.timeIndex() << endl;
        Info << "Time = " << runTime.timeName() << endl;
        
        dimensionedScalar dt = runTime.deltaT();

        for (int intIter = 0; intIter < lusgsIters; intIter++)
        {
            twoFluidFlux.computeFlux();

            //volScalarField rho = alpha1*rho1 + alpha2*rho2;
            //volVectorField virtualVelocity = (alpha1*rho1*U1 + alpha2*rho2*U2)/rho;

            //virtualMassTerm = 0.5*rho1*(1.0 - alpha2)*alpha2*(DU2 - DU1);

            volScalarField rezAlphaRho1(-dt*(
                fvc::ddt(conservative.alphaRho1()) + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg())
            ));

            volScalarField rezAlphaRho2(-dt*(
                fvc::ddt(conservative.alphaRho2()) + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg())
            ));

            volVectorField rezAlphaRhoU1(-dt*(
                fvc::ddt(conservative.alphaRhoU1())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha1_pos()*mesh.Sf(), twoFluidFlux.alpha1_neg()*mesh.Sf())
                + drag.K()*(fluid.U1() - fluid.U2())
                //+ virtualMassCoeffsMinus1
                //+ virtualMassTerm
            ));

            volVectorField rezAlphaRhoU2(-dt*(
                fvc::ddt(conservative.alphaRhoU2())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
                - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha2_pos()*mesh.Sf(), twoFluidFlux.alpha2_neg()*mesh.Sf())
                + drag.K()*(fluid.U2() - fluid.U1())
                //+ virtualMassCoeffsMinus2
                //- virtualMassTerm
            ));

            volScalarField rezEpsilon1(-dt*(
                fvc::ddt(conservative.epsilon1())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg())
                + ((drag.K()*(fluid.U1() - fluid.U2())) & U1)
                + (virtualMassCoeffsMinus1 & U1)
                //+ (virtualMassTerm & U1)
                //+ ((dragTerm /*+ virtualMassTerm*/) & virtualVelocity)
            ));

            volScalarField rezEpsilon2(-dt*(
                fvc::ddt(conservative.epsilon2())
                + TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg())
                + ((drag.K()*(fluid.U2() - fluid.U1())) & U2)
                //+ (virtualMassCoeffsMinus2 & U2)
                //- (virtualMassTerm & U2)
                //- ((dragTerm /*+ virtualMassTerm*/) & virtualVelocity)
            ));

            dp     = dimensionedScalar("pzero", dimPressure, 0.0);
            dalpha = dimensionedScalar("alphazero", dimless, 0.0);
            dpsi   = dimensionedScalar("psizero", dimless, 0.0);
            dU1    = dimensionedVector("U1zero", dimVelocity, vector(0, 0, 0));
            dU2    = dimensionedVector("U2zero", dimVelocity, vector(0, 0, 0));
            dT1    = dimensionedScalar("T1zero", dimTemperature, 0.0);
            dT2    = dimensionedScalar("T2zero", dimTemperature, 0.0);

            //Info << "ok" << endl;
            psi = log(alpha/(1.0 - alpha));
            //Info << "ok2" << endl;

            //Info << "ok" << endl; 

            #include "lusgsSweep.H"

            //Info << "ok5" << endl;

            p     += dp;
            psi   += dpsi;
            //alpha += dalpha;
            U1    += dU1;
            U2    += dU2;
            T1    += dT1;
            T2    += dT2;

            Info << "psi max:   " << max(psi) << "psi min:   " << min(psi) << endl;

            psi = min(max(psi, -36.0), 36.0);
            alpha = 1.0/(1.0 + exp(-psi));

            //Info << "alpha max: " << max(alpha) << "alpha min: " << min(alpha) << endl;

            fluid.blendVanishingFluid();
            //fluid.boundAlpha();
            //fluid.bound();
            fluid.correctBoundaryCondition();
            drag.correct();

            /*if (runTime.timeIndex() > 12570)
            {
                Info << ">>> Forcing write , alpha: <<<" << alpha[15248] << endl;
                runTime.writeNow();
                runTime.write();
            }*/

            fluid.correctThermo();
            rho1.ref() = thermo1.rho();
            rho2.ref() = thermo2.rho();

            if(intIter == lusgsIters - 1) fluid.correctInterfacialPressure();
            fluid.correctConservative();

            /*DU1 = fvc::ddt(U1) + fvc::div(fvc::flux(U1), U1);
            DU2 = fvc::ddt(U2) + fvc::div(fvc::flux(U2), U2);
            virtualMassCoeffs = 0.*0.5*thermo1.rho()*(1.0 - alpha2)*alpha2;
            virtualMassCoeffsMinus1 = virtualMassCoeffs*(fvc::ddt(U2) + fvc::div(fvc::flux(U2), U2) - fvc::div(fvc::flux(U1), U1));
            virtualMassCoeffsMinus2 = virtualMassCoeffs*(fvc::ddt(U1) + fvc::div(fvc::flux(U1), U1) - fvc::div(fvc::flux(U2), U2));*/

            scalar finalRezp     = fvc::domainIntegrate(mag(dp)    /dt).value();
            //scalar finalRezalpha = fvc::domainIntegrate(mag(dalpha)/dt).value();
            scalar finalRezpsi   = fvc::domainIntegrate(mag(dpsi)  /dt).value();
            scalar finalRezU1    = fvc::domainIntegrate(mag(dU1)   /dt).value();
            scalar finalRezU2    = fvc::domainIntegrate(mag(dU2)   /dt).value();
            scalar finalRezT1    = fvc::domainIntegrate(mag(dT1)   /dt).value();
            scalar finalRezT2    = fvc::domainIntegrate(mag(dT2)   /dt).value();

            Info << "LUSGS:  Solving for p,     " << "Final residual = " << finalRezp     << ", No Iterations = " << intIter << nl;
            //Info << "LUSGS:  Solving for alpha, " << "Final residual = " << finalRezalpha << ", No Iterations = " << intIter << nl;
            Info << "LUSGS:  Solving for psi,   " << "Final residual = " << finalRezpsi   << ", No Iterations = " << intIter << nl;
            Info << "LUSGS:  Solving for U1,    " << "Final residual = " << finalRezU1    << ", No Iterations = " << intIter << nl;
            Info << "LUSGS:  Solving for U2,    " << "Final residual = " << finalRezU2    << ", No Iterations = " << intIter << nl;
            Info << "LUSGS:  Solving for T1,    " << "Final residual = " << finalRezT1    << ", No Iterations = " << intIter << nl;
            Info << "LUSGS:  Solving for T2,    " << "Final residual = " << finalRezT2    << ", No Iterations = " << intIter << nl;

            /*bool lastIteration = ((intIter + 1) == lusgsIntIters);

            lastIteration = lastIteration || ( 
                (finalRezp     < lusgsTolerance) && 
                (finalRezalpha < lusgsTolerance) && 
                (finalRezU1    < lusgsTolerance) && 
                (finalRezU2    < lusgsTolerance) && 
                (finalRezT1    < lusgsTolerance) && 
                (finalRezT2    < lusgsTolerance) );

            if (lastIteration)
            {
                Info << "LUSGS: converged in " << intIter + 1 << " iterations." << nl << nl;
                break;
            }*/
        }



        runTime.write();

        /*if (runTime.timeIndex() > 6930)
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
