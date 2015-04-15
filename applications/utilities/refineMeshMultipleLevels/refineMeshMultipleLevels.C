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

Description
    Utility to refine cells in multiple directions.

    Either supply -all option to refine all cells (3D refinement for 3D
    cases; 2D for 2D cases) or reads a refineMeshDict with
    - cellSet to refine
    - directions to refine

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "undoableMeshCutter.H"
#include "hexCellLooper.H"
#include "cellSet.H"
#include "twoDPointCorrector.H"
#include "directions.H"
#include "OFstream.H"
#include "multiDirRefinement.H"
#include "labelIOList.H"
#include "wedgePolyPatch.H"
#include "plane.H"
#include "topoSetSource.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Max cos angle for edges to be considered aligned with axis.
static const scalar edgeTol = 1E-3;


// Calculate some edge statistics on mesh.
void printEdgeStats(const primitiveMesh& mesh)
{
    label nX = 0;
    label nY = 0;
    label nZ = 0;

    scalar minX = GREAT;
    scalar maxX = -GREAT;
    vector x(1, 0, 0);

    scalar minY = GREAT;
    scalar maxY = -GREAT;
    vector y(0, 1, 0);

    scalar minZ = GREAT;
    scalar maxZ = -GREAT;
    vector z(0, 0, 1);

    scalar minOther = GREAT;
    scalar maxOther = -GREAT;

    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        vector eVec(e.vec(mesh.points()));

        scalar eMag = mag(eVec);

        eVec /= eMag;

        if (mag(eVec & x) > 1-edgeTol)
        {
            minX = min(minX, eMag);
            maxX = max(maxX, eMag);
            nX++;
        }
        else if (mag(eVec & y) > 1-edgeTol)
        {
            minY = min(minY, eMag);
            maxY = max(maxY, eMag);
            nY++;
        }
        else if (mag(eVec & z) > 1-edgeTol)
        {
            minZ = min(minZ, eMag);
            maxZ = max(maxZ, eMag);
            nZ++;
        }
        else
        {
            minOther = min(minOther, eMag);
            maxOther = max(maxOther, eMag);
        }
    }

    Pout<< "Mesh edge statistics:" << endl
        << "    x aligned :  number:" << nX << "\tminLen:" << minX
        << "\tmaxLen:" << maxX << endl
        << "    y aligned :  number:" << nY << "\tminLen:" << minY
        << "\tmaxLen:" << maxY << endl
        << "    z aligned :  number:" << nZ << "\tminLen:" << minZ
        << "\tmaxLen:" << maxZ << endl
        << "    other     :  number:" << mesh.nEdges() - nX - nY - nZ
        << "\tminLen:" << minOther
        << "\tmaxLen:" << maxOther << endl << endl;
}


// Return index of coordinate axis.
label axis(const vector& normal)
{
    label axisIndex = -1;

    if (mag(normal & point(1, 0, 0)) > (1-edgeTol))
    {
        axisIndex = 0;
    }
    else if (mag(normal & point(0, 1, 0)) > (1-edgeTol))
    {
        axisIndex = 1;
    }
    else if (mag(normal & point(0, 0, 1)) > (1-edgeTol))
    {
        axisIndex = 2;
    }

    return axisIndex;
}


//- Returns -1 or cartesian coordinate component (0=x, 1=y, 2=z) of normal
//  in case of 2D mesh
label twoDNess(const polyMesh& mesh)
{
    const pointField& ctrs = mesh.cellCentres();

    if (ctrs.size() < 2)
    {
        return -1;
    }


    //
    // 1. All cell centres on single plane aligned with x, y or z
    //

    // Determine 3 points to base plane on.
    vector vec10 = ctrs[1] - ctrs[0];
    vec10 /= mag(vec10);

    label otherCellI = -1;

    for (label cellI = 2; cellI < ctrs.size(); cellI++)
    {
        vector vec(ctrs[cellI] - ctrs[0]);
        vec /= mag(vec);

        if (mag(vec & vec10) < 0.9)
        {
            // ctrs[cellI] not in line with n
            otherCellI = cellI;

            break;
        }
    }

    if (otherCellI == -1)
    {
        // Cannot find cell to make decent angle with cell0-cell1 vector.
        // Note: what to do here? All cells (almost) in one line. Maybe 1D case?
        return -1;
    }

    plane cellPlane(ctrs[0], ctrs[1], ctrs[otherCellI]);


    forAll(ctrs, cellI)
    {
        const labelList& cEdges = mesh.cellEdges()[cellI];

        scalar minLen = GREAT;

        forAll(cEdges, i)
        {
            minLen = min(minLen, mesh.edges()[cEdges[i]].mag(mesh.points()));
        }

        if (cellPlane.distance(ctrs[cellI]) > 1E-6*minLen)
        {
            // Centres not in plane
            return  -1;
        }
    }

    label axisIndex = axis(cellPlane.normal());

    if (axisIndex == -1)
    {
        return axisIndex;
    }


    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    //
    // 2. No edges without points on boundary
    //

    // Mark boundary points
    boolList boundaryPoint(mesh.allPoints().size(), false);

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        forAll(patch, patchFaceI)
        {
            const face& f = patch[patchFaceI];

            forAll(f, fp)
            {
                boundaryPoint[f[fp]] = true;
            }
        }
    }


    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (!boundaryPoint[e.start()] && !boundaryPoint[e.end()])
        {
            // Edge has no point on boundary.
            return -1;
        }
    }


    // 3. For all non-wedge patches: all faces either perp or aligned with
    //    cell-plane normal. (wedge patches already checked upon construction)

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        if (!isA<wedgePolyPatch>(patch))
        {
            const vectorField& n = patch.faceAreas();

            scalarField cosAngle = mag(n/mag(n) & cellPlane.normal());

            if (mag(min(cosAngle) - max(cosAngle)) > 1E-6)
            {
                // cosAngle should be either ~1 over all faces (2D front and
                // back) or ~0 (all other patches perp to 2D)
                return -1;
            }
        }
    }

    return axisIndex;
}

// Copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const topoSet& fromSet,
    const word& toName
)
{
    Info<< "Backing up " << fromName << " into " << toName << endl;

    topoSet backupSet(mesh, toName, fromSet);

    backupSet.write();
}


// Read and copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const word& toName
)
{
    topoSet fromSet(mesh, fromName, IOobject::READ_IF_PRESENT);

    backup(mesh, fromName, fromSet, toName);
}


// make cell set (original cellSet.C
void makeCellSet(word setName, word actionName, const dictionary& dict,const polyMesh& mesh)
{
	topoSetSource::setAction action = topoSetSource::toAction(actionName);
	
	word dictName(setName+"Dict");
	
    // Create topoSetSources
    PtrList<topoSetSource> topoSetSources
    (
        dict.lookup(dictName),
        topoSetSource::iNew(mesh)
    );


    // Load set to work
    autoPtr<topoSet> currentSetPtr(NULL);
    IOobject::readOption r;

    if ((action == topoSetSource::NEW) || (action == topoSetSource::CLEAR))
    {
        r = IOobject::NO_READ;

        backup(mesh, setName, setName + "_old");

        currentSetPtr.reset
        (
            new cellSet
            (
                mesh,
                setName,
                mesh.nCells()/10+1  // Reasonable size estimate.
            )
        );
    }
    else
    {
        r = IOobject::MUST_READ;

        currentSetPtr.reset
        (
            new cellSet
            (
                mesh,
                setName,
                r
            )
        );
    }

    topoSet& currentSet = currentSetPtr();

    Info<< "Set:" << currentSet.name()
        << "  Size:" << currentSet.size()
        << "  Action:" << actionName
        << endl;

    if ((r == IOobject::MUST_READ) && (action != topoSetSource::LIST))
    {
        // currentSet has been read so can make copy.
        backup(mesh, setName, currentSet, setName + "_old");
    }

    if (action == topoSetSource::CLEAR)
    {
        // Already handled above by not reading
    }
    else if (action == topoSetSource::INVERT)
    {
        currentSet.invert(currentSet.maxSize(mesh));
    }
    else if (action == topoSetSource::LIST)
    {
        currentSet.writeDebug(Info, mesh, 100);
        Info<< endl;
    }
    else if (action == topoSetSource::SUBSET)
    {
        // Apply topoSetSources to it to handle new/add/delete
        forAll(topoSetSources, topoSetSourceI)
        {
            // Backup current set.
            topoSet oldSet(mesh, currentSet.name() + "_old2", currentSet);

            currentSet.clear();

            topoSetSources[topoSetSourceI].applyToSet
            (
                topoSetSource::NEW,
                currentSet
            );

            // Combine new value of currentSet with old one.
            currentSet.subset(oldSet);
        }
    }
    else
    {
        // Apply topoSetSources to it to handle new/add/delete
        forAll(topoSetSources, topoSetSourceI)
        {
            topoSetSources[topoSetSourceI].applyToSet(action, currentSet);
        }
    }


    if (action != topoSetSource::LIST)
    {
        // Set has changed.

        // Sync across coupled patches.
        currentSet.sync(mesh);

        Info<< "Writing " << currentSet.name()
            << " (size " << currentSet.size() << ") to "
            << currentSet.instance()/currentSet.local()
               /currentSet.name()
            << endl << endl;

        currentSet.write();
    }
}

// refine cell set (taken from refineMesh)
void refineCellSet(word setName, const dictionary& refineDict,polyMesh& mesh, Time & runTime,const bool& overwrite)
{
	//Save oldinstance name for ovewriting
	const word oldInstance = mesh.pointsInstance();

    // List of cells to refine
	cellSet cells(mesh, setName);

    Pout<< "Read " << cells.size() << " cells from cellSet "
        << cells.instance()/cells.local()/cells.name()
        << endl << endl;

    labelList refCells = cells.toc();

    string oldTimeName(runTime.timeName());

    if (!overwrite)
    {
        runTime++;
    }


    // Multi-directional refinement (does multiple iterations)
    multiDirRefinement multiRef(mesh, refCells, refineDict);


    // Write resulting mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();


    // Get list of cell splits.
    // (is for every cell in old mesh the cells they have been split into)
    const labelListList& oldToNew = multiRef.addedCells();


    // Create cellSet with added cells for easy inspection
    /*cellSet newCells(mesh, "refinedCells", refCells.size());

    forAll(oldToNew, oldCellI)
    {
        const labelList& added = oldToNew[oldCellI];

        forAll(added, i)
        {
            newCells.insert(added[i]);
        }
    }

    Pout<< "Writing refined cells (" << newCells.size() << ") to cellSet "
        << newCells.instance()/newCells.local()/newCells.name()
        << endl << endl;

    newCells.write();*/




    //
    // Invert cell split to construct map from new to old
    //

    labelIOList newToOld
    (
        IOobject
        (
            "cellMap",
            runTime.timeName(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.nCells()
    );
    newToOld.note() =
        "From cells in mesh at "
      + runTime.timeName()
      + " to cells in mesh at "
      + oldTimeName;


    forAll(oldToNew, oldCellI)
    {
        const labelList& added = oldToNew[oldCellI];

        if (added.size())
        {
            forAll(added, i)
            {
                newToOld[added[i]] = oldCellI;
            }
        }
        else
        {
            // Unrefined cell
            newToOld[oldCellI] = oldCellI;
        }
    }

    Info<< "Writing map from new to old cell to "
        << newToOld.objectPath() << nl << endl;

    newToOld.write();


    // Some statistics.
    printEdgeStats(mesh);
}


// Main program:

int main(int argc, char *argv[])
{
    Foam::argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"

    printEdgeStats(mesh);

    bool overwrite = args.optionFound("overwrite");

    // Dictionary to control refinement
    Info<< "Refining according to refineMeshDict" << nl << endl;

	dictionary refineDict =
		IOdictionary
		(
		    IOobject
		    (
		        "refineMeshDict",
		        runTime.system(),
		        mesh,
		        IOobject::MUST_READ,
		        IOobject::NO_WRITE
		    )
		);
	
	// Dictionary of refinement sets
	Info<< "Reading MyCellSetsDict\n" << endl;
    IOdictionary MyCellSetsDict
    (
        IOobject
        (
            "MyCellSetsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    wordList setNames(MyCellSetsDict.lookup("names"));
    wordList actionNames(MyCellSetsDict.lookup("actions"));

    forAll(setNames,iSet)
    {
    	makeCellSet(setNames[iSet],actionNames[iSet],MyCellSetsDict,mesh);
    	mesh.readUpdate();
    	refineCellSet(setNames[iSet],refineDict,mesh,runTime,overwrite);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
