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

Info << nl << "Setting up RBF radial expansion matrices " << endl;
scalar outerRadius = 1.0;
scalar innerRadius = 0.1;
vector origin = vector::zero;
const scalar& pi = mathematicalConstant::pi;

//Read things if present
scalar motionAmplitude = dynamicMeshDict.lookupOrDefault("motionAmplitude",0.05);
scalar expansionFrequency = dynamicMeshDict.lookupOrDefault("motionFrequency",5.0);
scalar rotationFrequency = dynamicMeshDict.lookupOrDefault("rotationFrequency",expansionFrequency);
scalar rotationAmplitude = dynamicMeshDict.lookupOrDefault("rotationAmplitude",360.0)/180.0*pi;
bool rotation = readBool(dynamicMeshDict.lookup("rotation"));
bool expansion = readBool(dynamicMeshDict.lookup("expansion"));

if(rotation && expansion){
    Info << "Doing rotation AND expansion as motion." << endl;
}else if(rotation){
    Info << "Doing rotation as motion." << endl;
}else if(expansion){
    Info << "Doing expansion as motion." << endl;
}else{
    Info << "Doing no motion." << endl;
}

//Calculate things
scalarField radius = Foam::sqrt(Foam::sqr(mesh.points().component(0)()-origin.x()) + Foam::sqr(mesh.points().component(1)()-origin.y()));
scalarField theta = Foam::atan2(mesh.points().component(1)(),mesh.points().component(0)());

//Set all patches to fixed
labelList fixedPatches(mesh.boundaryMesh().size(),0);
label counter=0;
forAll(mesh.boundaryMesh(),iPatch){
    if(mesh.boundaryMesh()[iPatch].type()!="empty"){
        fixedPatches[counter]=iPatch;
        counter++;
    }
}
fixedPatches.setSize(counter);

List<scalarField> fixedPatchesRadius(fixedPatches.size());
forAll(fixedPatchesRadius,ipatch){
    const label& patchID = fixedPatches[ipatch];
    const labelList& meshP = mesh.boundaryMesh()[patchID].meshPoints();
    fixedPatchesRadius[ipatch] = scalarField(meshP.size(),0.0);
    forAll(meshP,ipoint){
        fixedPatchesRadius[ipatch][ipoint] = radius[meshP[ipoint]];
    }
}

//Control points
scalarField rc(3,0.0);
rc[0]=innerRadius;
rc[1]=innerRadius+(outerRadius-innerRadius)*0.5;
rc[2]=outerRadius;

scalarSquareMatrix controlA(3,0.0);
for(int i=0;i<3;i++){
    for(int j=i+1;j<3;j++){
        scalar r=mag(rc[i]-rc[j]);
        if(r>SMALL){
            controlA[i][j]=Foam::sqr(r)*Foam::log(r);
            controlA[j][i]=Foam::sqr(r)*Foam::log(r);
        }
    }
}

LUscalarMatrix solveControlA(controlA);
Info << "Build RBF LU Matrix" << endl;

RectangularMatrix<scalar> evaluationH(mesh.points().size(),3,0.0);
for(int i=0;i<mesh.points().size();i++){
    for(int j=0;j<3;j++){
        scalar r=mag(radius[i]-rc[j]);
        if(r>SMALL){
            evaluationH[i][j]=Foam::sqr(r)*Foam::log(r);
        }
    }
}
Info << "Matrix H build with size = " << mesh.points().size()*3*8/1000/1000 << "MB" << nl << endl;

scalarField motionRc(3,0.0);
