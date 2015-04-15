vectorField movingControlPointsMotion = vectorField(movingControlPoints.size(),vector::zero);
const scalar& t = runTime.time().value();
const scalar& dt = runTime.deltaT().value();

//Set smoothing factor
scalar smoothingFactor = 1;
scalar smoothingFactorPrev = 1;
if(t < 1.0/motionFreq)
{
    smoothingFactor = 0.5 - 0.5*Foam::cos(pi*motionFreq*t);
}
if((t-dt) < 1.0/motionFreq)
{
    smoothingFactorPrev = 0.5 - 0.5*Foam::cos(pi*motionFreq*(t-dt));
}
//Info << "Smoothing new|prev = " << smoothingFactor << "|" << smoothingFactorPrev << endl;

//Set translation
vector translation = motionAmplitude*(smoothingFactor*Foam::sin(2*pi*motionFreq*t)-smoothingFactorPrev*Foam::sin(2*pi*motionFreq*(t-dt)));
//Info << "translation = " << translation << endl;

//Set rotation
scalar rotAngle = rotationAmplitude*smoothingFactor*Foam::sin(2*pi*motionFreq*t);
scalar rotAnglePrev = rotationAmplitude*smoothingFactorPrev*Foam::sin(2*pi*motionFreq*(t-dt));
tensor RzCur(Foam::cos(rotAngle), -Foam::sin(rotAngle), 0, Foam::sin(rotAngle), Foam::cos(rotAngle), 0, 0, 0, 1);
tensor RzCurPrev(Foam::cos(rotAnglePrev), -Foam::sin(rotAnglePrev), 0, Foam::sin(rotAnglePrev), Foam::cos(rotAnglePrev), 0, 0, 0, 1);
//Info << "rotAngle new|prev = " << rotAngle << "|" << rotAnglePrev << endl;

vectorField rotation = (RzCur & ( movingControlPoints - rotationOrigin )) + rotationOrigin;
vectorField rotationPrev = (RzCurPrev & ( movingControlPoints - rotationOrigin )) + rotationOrigin;

//Gather deformation
movingControlPointsMotion = translation + (rotation - rotationPrev);

vectorField newPoints = rbfIntRed.interpolateRed(movingControlPointsMotion);
newPoints += mesh.points();
mesh.movePoints(newPoints);
//bool meshChanged = mesh.update();
bool meshChanged = true;
