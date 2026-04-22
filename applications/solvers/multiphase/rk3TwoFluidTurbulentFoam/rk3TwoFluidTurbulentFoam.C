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
#include "PhaseCompressibleTurbulenceModel.H"
#include "ThermalDiffusivity.H"
#include "fvcSmooth.H"

#include "processorPolyPatch.H"

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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    turbulence1->correct();

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        amaxSf = fluid.amaxSf();
        #include "courantNo.H"
        #include "setDeltaT.H"
        ++runTime;

        Info << nl << "Iter = " << runTime.timeIndex() << endl;
        Info<< "Time = " << runTime.timeName() << endl;

        dimensionedScalar dt = runTime.deltaT();

        scalar coeff[3][3] = {{1.0, 0.0 , 1.0}, {3.0/4.0, 1.0/4.0, 1.0/4.0}, {1.0/3.0, 2.0/3.0, 2.0/3.0}};

        for (size_t i = 0; i < 3; i++)
        {
            twoFluidFlux.computeFlux();

            volScalarField dragK = drag.K();

            conservative.alphaRho1() = coeff[i][0]*conservative.alphaRho1().oldTime()
                                     + coeff[i][1]*conservative.alphaRho1()
                                     - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux1_pos(), twoFluidFlux.alphaRhoFlux1_neg()));
            
            conservative.alphaRho2() = coeff[i][0]*conservative.alphaRho2().oldTime()
                                     + coeff[i][1]*conservative.alphaRho2()
                                     - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoFlux2_pos(), twoFluidFlux.alphaRhoFlux2_neg()));

            conservative.alphaRhoU1() = coeff[i][0]*conservative.alphaRhoU1().oldTime()
                                      + coeff[i][1]*conservative.alphaRhoU1()
                                      - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux1_pos(), twoFluidFlux.alphaRhoUFlux1_neg())
                                          + alpha1*fvc::div(turbulence1->devRhoReff())
                                          - fluid.pInt()*TwoFluidFoam::fvc::div((1.0 - twoFluidFlux.alpha2_pos())*mesh.Sf(), (1.0 - twoFluidFlux.alpha2_neg())*mesh.Sf())
                                          + dragK*(fluid.U1() - fluid.U2()));
            
            conservative.alphaRhoU2() = coeff[i][0]*conservative.alphaRhoU2().oldTime()
                                      + coeff[i][1]*conservative.alphaRhoU2()
                                      - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoUFlux2_pos(), twoFluidFlux.alphaRhoUFlux2_neg())
                                          - fluid.pInt()*TwoFluidFoam::fvc::div(twoFluidFlux.alpha2_pos()*mesh.Sf(), twoFluidFlux.alpha2_neg()*mesh.Sf())
                                          + dragK*(fluid.U2() - fluid.U1()));

            conservative.epsilon1() = coeff[i][0]*conservative.epsilon1().oldTime()
                                    + coeff[i][1]*conservative.epsilon1()
                                    - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux1_pos(), twoFluidFlux.alphaRhoEFlux1_neg())
                                    + alpha1*fvc::div(turbulence1->devRhoReff() & U1)
                                    - alpha1*fvc::laplacian(turbulence1->alphaEff(), h1)
                                    //- fvc::laplacian(turbulence1->kappaEff(), T1)
                                    + ((dragK*(fluid.U1() - fluid.U2())) & U1));
            
            conservative.epsilon2() = coeff[i][0]*conservative.epsilon2().oldTime()
                                    + coeff[i][1]*conservative.epsilon2()
                                    - coeff[i][2]*dt*(TwoFluidFoam::fvc::div(twoFluidFlux.alphaRhoEFlux2_pos(), twoFluidFlux.alphaRhoEFlux2_neg())
                                    + ((dragK*(fluid.U2() - fluid.U1())) & U2));

            fluid.correct();
            fluid.blendVanishingFluid();
            fluid.correctBoundaryCondition();
            fluid.correctThermo();
            if(i == 2) fluid.correctInterfacialPressure();
            fluid.correctConservative();

            drag.correct();
        }

        scalar finalRezp     = fvc::domainIntegrate(mag(p - p.oldTime())        /dt).value();
        scalar finalRezalpha = fvc::domainIntegrate(mag(alpha2 - alpha2.oldTime())/dt).value();
        scalar finalRezU1    = fvc::domainIntegrate(mag(U1 - U1.oldTime())      /dt).value();
        scalar finalRezU2    = fvc::domainIntegrate(mag(U2 - U2.oldTime())      /dt).value();
        scalar finalRezT1    = fvc::domainIntegrate(mag(T1 - T1.oldTime())      /dt).value();
        scalar finalRezT2    = fvc::domainIntegrate(mag(T2 - T2.oldTime())      /dt).value();

        Info << "rk3:  Solving for p,     " << "Final residual = " << finalRezp     << nl;
        Info << "rk3:  Solving for alpha, " << "Final residual = " << finalRezalpha << nl;
        Info << "rk3:  Solving for U1,    " << "Final residual = " << finalRezU1    << nl;
        Info << "rk3:  Solving for U2,    " << "Final residual = " << finalRezU2    << nl;
        Info << "rk3:  Solving for T1,    " << "Final residual = " << finalRezT1    << nl;
        Info << "rk3:  Solving for T2,    " << "Final residual = " << finalRezT2    << nl;

        turbulence1->correct();
        h1 = thermo1.he() + p/thermo1.rho();
        phi1 = linearInterpolate(U1) & mesh.Sf();
        U = U1;

        rho1.ref() = thermo1.rho();
        rho2.ref() = thermo2.rho();

        runTime.write();

        #include "../tools/boundaryFlux.H"

        //runTime.printExecutionTime(Info);

        //std::cin.ignore();
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
