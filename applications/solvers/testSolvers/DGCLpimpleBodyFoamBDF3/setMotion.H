//Get motion of bodies and update mesh accordingly
//bodyCol.update();
//ms.setMotion(bodyCol.getPatchMotion());
// Call interpolation
const scalar& pi = mathematicalConstant::pi;
vectorField movingControlPointsMotion = vectorField(movingControlPoints.size(),vector::zero);
movingControlPointsMotion = motionVector*(Foam::sin(2*pi*motionFrequency*runTime.value())-Foam::sin(2*pi*motionFrequency*(runTime.value()-runTime.deltaT().value())));

vectorField newPoints = rbfIntRed.interpolateRed(movingControlPointsMotion);
newPoints += mesh.points();
mesh.movePoints(newPoints);
//bool meshChanged = mesh.update();
bool meshChanged = true;
