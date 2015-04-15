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
    argList::validArgs.append("radius");
    argList::validOptions.insert("overwrite", "");

    #   include "addTimeOptions.H"

    #   include "setRootCase.H"
    #   include "createTime.H"
    runTime.functionObjects().off();

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
    #   include "checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);

    #   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();
    bool overwrite = args.optionFound("overwrite");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //Get patches
    wordList patches(IStringStream(args.additionalArgs()[0])());
    scalarField radius(IStringStream(args.additionalArgs()[1])());
    Info << "patches = " << patches << endl;
    Info << "radius = " << radius << endl;

    if(patches.size()!=radius.size()){
        FatalErrorIn("main writePatchCoordinates")
        << "patches hase size " << patches.size() << "!= size of radius = " << radius.size()
        << abort(FatalError);
    }

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

    vectorField newpoints = mesh.points();//create newpoints
    forAll(patchIDs,ipatch){
        const label& patchID = patchIDs[ipatch];
        const scalar& desRadius = radius[ipatch];
        const labelList& meshP = mesh.boundaryMesh()[patchID].meshPoints();

        scalarField patchRadius(meshP.size(),0.0);
        forAll(patchRadius,iPoint){
            const vector& patchPoint = mesh.points()[meshP[iPoint]];
            patchRadius[iPoint] = Foam::sqrt(Foam::sqr(patchPoint.x()) + Foam::sqr(patchPoint.y()));
        }
        Info << "maxError in outer radius on patch " << patchID << " = " << max(mag((patchRadius - desRadius)/desRadius)) << endl;

        //Correct radial position patch
        scalarField thetas(meshP.size(),0.0);
        forAll(patchRadius,iPoint){
            const vector& patchPoint = mesh.points()[meshP[iPoint]];
            thetas[iPoint]=Foam::atan2(patchPoint.y(),patchPoint.x());
            newpoints[meshP[iPoint]].x()=desRadius*Foam::cos(thetas[iPoint]);
            newpoints[meshP[iPoint]].y()=desRadius*Foam::sin(thetas[iPoint]);
        }
    }
    mesh.movePoints(newpoints);

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    Info<< "Writing mesh to " << oldInstance << endl;
    mesh.write();

    forAll(patchIDs,ipatch){
        const label& patchID = patchIDs[ipatch];
        const scalar& desRadius = radius[ipatch];
        const labelList& meshP = mesh.boundaryMesh()[patchID].meshPoints();

        scalarField patchRadius(meshP.size(),0.0);
        forAll(patchRadius,iPoint){
            const vector& patchPoint = mesh.points()[meshP[iPoint]];
            patchRadius[iPoint] = Foam::sqrt(Foam::sqr(patchPoint.x()) + Foam::sqr(patchPoint.y()));
        }
        Info << "maxError in outer radius on patch " << patchID << " after correction = " << max(mag((patchRadius - desRadius)/desRadius)) << endl;
    }
}


// ************************************************************************* //
