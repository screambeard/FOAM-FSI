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
    Transient solver for incompressible, laminar flow of Newtonian fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "gust.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "bodyCollector.H"
#include "RBFMotionSolverExt.H"

#include "processorFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "readPIMPLEControls.H"
#   include "createFields.H"
#	include "readGustProperties.H"
#   include "initContinuityErrs.H"
#   include "readBodyProperties.H"//Reading body properties
//#	include "checkTimeDiscretisationSchemeRC.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "CourantNo.H"

		//Make absolute phi (no gust)
		phi = phi - phiGust;//gust
        fvc::makeAbsolute(phi, U);//mesh
	    
#       include "setDeltaT.H"

		runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Get motion of bodies and update mesh accordingly
        bodyCol.update();
        ms.setMotion(bodyCol.getPatchMotion());
        bool meshChanged = mesh.update();
#       include "updateSf.H"

		//Inlcude gust fluxes if required
#		include "setGust.H"
#		include "setArtificialVolumeChange.H"

        if (correctPhi && (mesh.moving() || meshChanged))
        {
#           include "correctPhi.H"
        }

        //Make relative
        fvc::makeRelative(phi,U);
        phi = phi + phiGust;

        if (mesh.moving() && checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

        //Check DGCL
#       include "checkDGCL.H"

		// --- Pressure-velocity PIMPLE corrector loop
		for (int ocorr=0; ocorr<nOuterCorr; ocorr++)
		{
			if (nOuterCorr != 1)
            {
                p.storePrevIter();
            } 

#			include "UEqn.H"

		    // --- PISO loop
		    for (int corr=0; corr<nCorr; corr++)
		    {
#				include "pEqn.H"
		    }

            //Make relative phi
		    turbulence->correct();
		}
		
		Utotal == U + fvc::average(Ugust);
		Utotal.rename("Utotal");
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
