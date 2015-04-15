/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    pimpleDyMFoam.C

Description
    Transient solver for incompressible, flow of Newtonian fluids
    on a moving mesh using the PIMPLE (merged PISO-SIMPLE) algorithm.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"
#include "bodyCollector.H"
#include "RBFMotionSolverExt.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "readBodyProperties.H"//Reading body properties

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        fvc::makeAbsolute(phi, U);

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Get motion of bodies and update mesh accordingly
        bodyCol.update();
        ms.setMotion(bodyCol.getPatchMotion());
        bool meshChanged = mesh.update();
#       include "updateSf.H"

        if (correctPhi && (mesh.moving() || meshChanged))
        {
#           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);
        if (mesh.moving() && checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

        label oCorr=1;
        scalar relativeResidual = 1;
        bool converged=false;
        bool lastIter=false;
        do
        {
            Info << "outer iteration: " << oCorr << endl;
            maxResidual = 0;//reset maxResidual for new outerloop

            p.storePrevIter();
            U.storePrevIter();
#           include "UEqn.H"

            // --- PISO loop
            label corr=0;
            do
            {
                corr++;
#               include "pEqn.H"
            }while(eqnResidual > innerConvergence && corr < nCorr);

#           include "movingMeshContinuityErrs.H"

            //Info << "Doing Relexation on turbulence" << endl;
            turbulence->correct();

            //Get residuals
            scalar residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / gSumMag( p.internalField() );
            scalar residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / gSumMag( U.internalField() );
            relativeResidual = max( residualPressure, residualVelocity );
            Info << "residualPressure = " << residualPressure << ", residualVelocity = " << residualVelocity << endl;

            //Check if this was the last iter and set converged to true
            if(lastIter){
                converged = true;
            }
            //Check convergence and set lastIter to true
            if(relativeResidual < outerConvergence || oCorr == nOuterCorr){
                lastIter = true;
            }

            oCorr++;
            //Info << "lastIter|converged = " << lastIter << "|" << converged << endl;
        }while(!converged);


        runTime.write();
        //Write bodyMotion properties to file (if requested)
		bodyCol.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
