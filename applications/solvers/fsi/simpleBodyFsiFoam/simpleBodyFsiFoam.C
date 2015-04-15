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
    icoFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids with
    mesh motion.  Set up as a fake fluid structure interaction solver

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"
    
#include "RBFMotionSolverExt.H"
#include "fsiInterface.H"
#include "bodyCollector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"//Set location of case
#   include "createTime.H"//Create runtime environment
#   include "createDynamicFvMesh.H"//Create dynamic fluid mesh
#   include "createFields.H"//Create fluid fields for incompressible flow & read (constant) density (rhoFluid)
#   include "readBodyProperties.H"//Reading body properties: spring stifnesses, mass properties, cg, rc
#   include "readCouplingProperties.H"//Read coupling patches solid and fluid. Construct interpolators (fluid-solid, solid-fluid), mesh motion object (tForce, motionUFluidPatch, tppi)

#   include "initContinuityErrs.H"//Initialize continuity error calculations
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        p.storePrevIter();

		//Update old values because new time step
        bodyCol.update();

		//Re-initialize variables each time step
		fsiInter.reset();

		do{
        	Info<< "Fsi iteration = " << fsiInter.iter() << endl;

	#       include "setMotion.H"
	#       include "solveFluid.H"
			
			fsiInter.update(rhoFluid,nu,U,p);
		} while(!fsiInter.converged());

		//Writing data
        runTime.write();
        //Write bodyMotion properties to file (if requested)
		bodyCol.write();
		//Write fsi properties to file (if requested)
		fsiInter.write();
		
        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }
    
    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
