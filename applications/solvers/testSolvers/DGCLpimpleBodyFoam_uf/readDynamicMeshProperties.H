IOdictionary dynamicMeshDict
(
    IOobject
    (
        "dynamicMeshDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
dynamicMeshDict.subDict("interpolation").set("Dim","ThreeD");

//Set motion variables
scalar motionFrequency = 10.0;
vector motionVector(0.2,0.2,0.0);

//Create all point markers for building matrices
labelList pointsMarker(mesh.points().size(),0);
vectorField controlPoints(0,vector::zero);
label counter=0;
forAll(mesh.boundaryMesh(),iPatch){
    if(mesh.boundaryMesh()[iPatch].type()!="empty"){
        const labelList& patchPointIDs = mesh.boundaryMesh()[iPatch].meshPoints();
        controlPoints.setSize(controlPoints.size()+patchPointIDs.size());
        forAll(patchPointIDs,ipoint){
            if(pointsMarker[patchPointIDs[ipoint]] == 0){
                controlPoints[counter] = mesh.points()[patchPointIDs[ipoint]];
                pointsMarker[patchPointIDs[ipoint]] = 1;
                counter++;
            }
        }
    }
}

//Set control points
controlPoints.setSize(counter+8);
controlPoints[counter] = vector(0.4,0.4,0);
controlPoints[counter+1] = vector(0.6,0.4,0);
controlPoints[counter+2] = vector(0.6,0.6,0);
controlPoints[counter+3] = vector(0.4,0.6,0);
controlPoints[counter+4] = vector(0.4,0.4,1);
controlPoints[counter+5] = vector(0.6,0.4,1);
controlPoints[counter+6] = vector(0.6,0.6,1);
controlPoints[counter+7] = vector(0.4,0.6,1);
vectorField intPoints = mesh.points();
Info << "controlPoints.size()  = " << controlPoints.size() << endl;

vectorField movingControlPoints(SubList<vector>(controlPoints,8,counter));
Info << "movingcontrolPoints = " << movingControlPoints << endl;

//Create RBF interpolation object
RBFInterpolationReduced rbfIntRed(mesh,dynamicMeshDict.subDict("interpolation"),controlPoints,intPoints);

//Create moving Control Point indices
labelList movingControlPointsIndex(8,0);
forAll(movingControlPointsIndex,i){
    movingControlPointsIndex[i] = counter+i;
}
rbfIntRed.createReducedEvaluationMatrix(movingControlPointsIndex);

//Create file for maximum p output
OFstream maxPoutput("maximumPressureInTime.dat");
maxPoutput << "t" << "\t" << "maxP" << endl;
