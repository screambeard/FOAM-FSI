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

forAll(thetas,iPoint)
{
    scalar theta=thetas[iPoint]/180.0*pi;

    vector deformation = vector::zero;
    deformation.x()=motionAmplitude*Foam::cos(theta*6)*Foam::cos(theta)*(smoothingFactor*Foam::sin(2*pi*motionFreq*t)-smoothingFactorPrev*Foam::sin(2*pi*motionFreq*(t-dt)));
    deformation.y()=motionAmplitude*Foam::cos(theta*6)*Foam::sin(theta)*(smoothingFactor*Foam::sin(2*pi*motionFreq*t)-smoothingFactorPrev*Foam::sin(2*pi*motionFreq*(t-dt)));

    movingControlPointsMotion[iPoint]=deformation;
    movingControlPointsMotion[iPoint+Nc]=deformation;
}

vectorField newPoints = rbfIntRed.interpolateRed(movingControlPointsMotion);
newPoints += mesh.points();
mesh.movePoints(newPoints);
bool meshChanged = true;
