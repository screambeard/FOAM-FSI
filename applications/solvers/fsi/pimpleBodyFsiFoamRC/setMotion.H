	//Get motion based on forces from fsiInterface. fsiInter.getPatchForces() has forces in order of movingPatchIDs_ from bodyCollector
    ms.setMotion(bodyCol.getPatchMotion(fsiInter.getPatchForces()));
    bool meshChanged = mesh.update();

#   include "updateSf.H"
#	include "checkVolContinuity.H"//Check moving mesh continuity errors