// vectorField movingControlPointsMotion = vectorField(movingControlPoints.size(),vector::zero);
Info << nl << "Setting displacment for radial expansion" << endl;
const scalar& t = runTime.time().value();
const scalar& dt = runTime.deltaT().value();
vectorField newPoints(mesh.points().size(),vector::zero);
newPoints = mesh.points();

// === Radial expansion calculation ==== //
if(expansion){
    scalar r = motionAmplitude*0.5*(1-Foam::cos(4*pi*t))*Foam::sin(2*pi*t);
    scalar r0 = motionAmplitude*0.5*(1-Foam::cos(4*pi*(t-dt)))*Foam::sin(2*pi*(t-dt));
    motionRc[1] = r-r0;
    //motionRc[1]=motionAmplitude*(Foam::sin(2*pi*expansionFrequency*t)-Foam::sin(2*pi*expansionFrequency*(t-dt)));

    scalarField beta=motionRc;
    Info << "Start solving for beta" << endl;
    solveControlA.solve(beta);
    Info << "Solved for beta. Beta = " << beta << endl;

    scalarField deltaRadius(mesh.points().size(),0.0);
    for(int i=0;i<mesh.points().size();i++){
        scalar dradius = 0.0;
        for(int j=0;j<3;j++){
            dradius += evaluationH[i][j]*beta[j];
        }
        radius[i] += dradius;
    }
}

// === Rotation calculation ==== //
if(rotation){
    //scalar dtheta = rotationAmplitude*rotationFrequency*dt;
    scalar dtheta = rotationAmplitude/2.0*(Foam::cos(pi*rotationFrequency*(t-dt))-Foam::cos(pi*rotationFrequency*t));
    forAll(theta,ipoint){
        theta[ipoint] += dtheta;
    }
}

//calculate new points
bool meshChanged = false;
if(expansion || rotation){
    forAll(radius,ipoint){
        newPoints[ipoint].x() = radius[ipoint]*Foam::cos(theta[ipoint]);
        newPoints[ipoint].y() = radius[ipoint]*Foam::sin(theta[ipoint]);
    }

    //Set all points of boundary to original position
    forAll(fixedPatches,iPatch){
        label patchID=fixedPatches[iPatch];
        const labelList& pointIDs = mesh.boundaryMesh()[patchID].meshPoints();
        forAll(pointIDs,ipoint){
            newPoints[pointIDs[ipoint]].x() = fixedPatchesRadius[iPatch][ipoint]*Foam::cos(theta[pointIDs[ipoint]]);
            newPoints[pointIDs[ipoint]].y() = fixedPatchesRadius[iPatch][ipoint]*Foam::sin(theta[pointIDs[ipoint]]);
        }
    }

    mesh.movePoints(newPoints);
    meshChanged = true;
}
