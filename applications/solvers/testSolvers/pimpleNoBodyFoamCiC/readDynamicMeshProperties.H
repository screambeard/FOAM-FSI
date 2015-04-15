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

//Set radial displacement control points
const scalar& pi = mathematicalConstant::pi;
label Nc=12;
scalar dTheta=360.0/Nc;
scalarField thetas(Nc,0.0);
vectorField xyMc(Nc,vector::zero);
forAll(thetas,itheta){
    thetas[itheta]=itheta*dTheta;
    xyMc[itheta].x()=0.5*Foam::cos(thetas[itheta]/180.0*pi);
    xyMc[itheta].y()=0.5*Foam::sin(thetas[itheta]/180.0*pi);
}

controlPoints.setSize(counter+2*Nc);
forAll(xyMc,iMc){
    controlPoints[counter+iMc]=xyMc[iMc];
    vector tmpVector=xyMc[iMc];
    tmpVector.z()=1.0;
    controlPoints[counter+iMc+Nc]=tmpVector;
}

vectorField intPoints = mesh.points();
Info << "controlPoints.size()  = " << controlPoints.size() << endl;

vectorField movingControlPoints(SubList<vector>(controlPoints,2*Nc,counter));
Info << "movingcontrolPoints = " << movingControlPoints << endl;

RBFInterpolationReduced rbfIntRed(mesh,dynamicMeshDict.subDict("interpolation"),controlPoints,intPoints);

labelList movingControlPointsIndex(2*Nc,0);
forAll(movingControlPointsIndex,i){
    movingControlPointsIndex[i] = counter+i;
}
rbfIntRed.createReducedEvaluationMatrix(movingControlPointsIndex);

scalar motionFreq = 1.0/0.2;
scalar motionAmplitude = 0.1;
