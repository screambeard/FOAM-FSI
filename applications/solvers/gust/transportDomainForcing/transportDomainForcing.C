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

#include "processorFvsPatchFields.H"
#include "inletOutletFvPatchFields.H"

#include "cellSet.H"
#include "faceSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#	include "readGustProperties.H"
#   include "initContinuityErrs.H"
#	include "checkTimeDiscretisationSchemeRC.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"
    
		runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

		//Calculate gust body force
#		include "setGust.H"

        if(correctPhi){
#			include "correctPhi.H"        	
        }

		for (int outCorr=0; outCorr<nOuterCorr; outCorr++)
		{

		    fvVectorMatrix UEqn
		    (
		        fvm::ddt(U)
		      + fvm::div(phi, U)
		    );

		    solve(UEqn == gustSource);
		}
		
		phi = ((fvc::interpolate(U)) & mesh.Sf());
		forAll(phi.boundaryField(),iPatch)
		{
			if(!phi.boundaryField()[iPatch].coupled())
			{
				phi.boundaryField()[iPatch] = U.boundaryField()[iPatch] & mesh.Sf().boundaryField()[iPatch];
			}
		}
		//Set phi at faces equal to left (upwind)
		forAll(forceFaceSet,iFace)
		{
			label faceID = forceFaceSet.toc()[iFace];
			phi.internalField()[faceID] = (U.internalField()[leftCellIDs[iFace]] & mesh.Sf().internalField()[faceID]);
		}

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
