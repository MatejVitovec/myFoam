/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: Open Source CFD
   \\    /   O peration     | 
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  | 
-------------------------------------------------------------------------------
License
    This file isn't part of foam-extend nor OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockLusgsFoam

Description
    Density-based compressible implicit steady-state & transient flow solver
    using LU-SGS method.

Author
    Jiri Furst
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "fvOptions.H"
#include "psiThermo.H"

#include "convectiveFlux.H"

#include <eigen3/Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " Riemann problem flux calculation and blocLUSGS implicit scheme"
    );

    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createTimeControls.H"

    #include "readLUSGSControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    while (runTime.run())
    {
        #include "readTimeControls.H"

        amaxSf = mag(fvc::interpolate(U) & mesh.Sf()) + mesh.magSf()*fvc::interpolate(a);
        #include "courantNo.H"
        #include "setDeltaT.H"
        //if (LTS)
        //{
        //    #include "setRDeltaTau.H"
        //}
        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl;
        Info<< "Iteration = " << runTime.timeIndex() << nl;
        
        scalar initialRezRho = 0, initialRezRhoU = 0, initialRezRhoE = 0;

        for (int intIter = 0; intIter < lusgsIntIters; intIter++)
        {
            Info << "LUSGS: iteration " << intIter + 1 << nl;
            
            flux.computeFlux();
            
            dimensionedScalar dt = runTime.deltaT();
            
            volScalarField dRho(
                -dt*(fvc::ddt(rho) + fvc::div(flux.rhoFlux()))
            );
            
            volVectorField dRhoU(
                -dt*(fvc::ddt(rhoU) + fvc::div(flux.rhoUFlux()))
            ); 

            volScalarField dRhoE(
                -dt*(fvc::ddt(rhoE) + fvc::div(flux.rhoEFlux()))
            );
            
            scalar rezRho  = fvc::domainIntegrate(mag(dRho) /dt).value();
            scalar rezRhoU = fvc::domainIntegrate(mag(dRhoU)/dt).value();
            scalar rezRhoE = fvc::domainIntegrate(mag(dRhoE)/dt).value();

            if (intIter == 0)
            {
                initialRezRho  = rezRho;
                initialRezRhoU = rezRhoU;
                initialRezRhoE = rezRhoE;
            }

            #include "blockLusgsSweep.H"

            rho  += dRho;
            rhoU += dRhoU;
            rhoE += dRhoE;
            
            #include "updateFields.H"

            scalar finalRezRho  = fvc::domainIntegrate(mag(dRho) /dt).value();
            scalar finalRezRhoU = fvc::domainIntegrate(mag(dRhoU)/dt).value();
            scalar finalRezRhoE = fvc::domainIntegrate(mag(dRhoE)/dt).value();

            Info << "LUSGS:  Solving for rho,  " 
                 << "Initial residual = " << rezRho << ", "
                 << "Final residual = " << finalRezRho << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for rhoU, " 
                 << "Initial residual = " << rezRhoU << ", "
                 << "Final residual = " << finalRezRhoU << ", No Iterations 1" << nl;
            Info << "LUSGS:  Solving for rhoE, " 
                 << "Initial residual = " << rezRhoE << ", "
                 << "Final residual = " << finalRezRhoE << ", No Iterations 1" << nl;

            bool lastIteration = ( intIter + 1 == lusgsIntIters );

            lastIteration = lastIteration || ( 
                (finalRezRho  < lusgsRelTol*initialRezRho) && 
                (finalRezRhoU < lusgsRelTol*initialRezRhoU) && 
                (finalRezRhoE < lusgsRelTol*initialRezRhoE));

            lastIteration = lastIteration || ( 
                (finalRezRho  < lusgsTolerance) && 
                (finalRezRhoU < lusgsTolerance) && 
                (finalRezRhoE < lusgsTolerance));

            if (lastIteration)
            {
                Info << "LUSGS: converged in " << intIter + 1 << " iterations." << nl;
                break;
            }
        }
	
        runTime.write();
	
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s\n" << endl;
    }

    Info << "\n end \n";

    return(0);
}


// ************************************************************************* //