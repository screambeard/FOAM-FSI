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
#include "fvMesh.H"  
#include "twoDPointCorrector.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("patches");
    timeSelector::addOptions();

#   include "setRootCase.H"//Set location of case
#   include "createTime.H"//Create runtime environment
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"//Create dynamic fluid mesh

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Determine whether 2D or not
    twoDPointCorrector twoDCorr(mesh);

    //Get patches
    wordList patches(IStringStream(args.additionalArgs()[0])());

    labelList patchIDs(patches.size(),-1);
    label nPoints = 0;
    forAll(patches,ipatch){
        label patchID = mesh.boundaryMesh().findPatchID(patches[ipatch]);
        if(patchID < 0){
            FatalErrorIn("main writePatchCoordinates")
            << "patchName " << patches[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs[ipatch] = patchID;
            nPoints += mesh.boundaryMesh()[patchID].localPoints().size();
        }
    }
    
    scalar avEmptyDir=0;
    if(twoDCorr.required()){
        nPoints /= 2;
        const vectorField& localpoints = mesh.boundaryMesh()[0].localPoints();
        scalar maxDir = max((twoDCorr.planeNormal() & localpoints));
        scalar minDir = min((twoDCorr.planeNormal() & localpoints));
        avEmptyDir = 0.5*(maxDir + minDir);
    }

    Info << "Total number of points for patches: " << nPoints << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        mesh.readUpdate();
        Info<< "\nTime = " << runTime.timeName() << endl;

        //Gather points
        vectorField patchPoints(nPoints,vector::zero);
        label counter=0;
        forAll(patchIDs,ipatch){
            label patchID = patchIDs[ipatch];
            vectorField localPoints = mesh.boundaryMesh()[patchID].localPoints();
            forAll(localPoints,ipoint){
                //Check whether in 2D only single side is picked
                if(twoDCorr.required()){
                    if((localPoints[ipoint] & twoDCorr.planeNormal()) > avEmptyDir){
                        patchPoints[counter] = localPoints[ipoint];
                        counter++;        
                    }
                }else{
                    patchPoints[counter] = localPoints[ipoint];
                    counter++;        
                }
            }
        }
        
        //Writing out in file
        OFstream of(fileName("patchCoordinates_" + timeDirs[timeI].name() + ".dat"));
        forAll(patchPoints,ipoint){
            if(ipoint==nPoints-1){
                of << patchPoints[ipoint].x() << " " << patchPoints[ipoint].y();
            }else{
                of << patchPoints[ipoint].x() << " " << patchPoints[ipoint].y() << endl;
            }
            
        }
    }
}


// ************************************************************************* //
