const scalar& time = runTime.time().value();

forAll(mesh.boundaryMesh(),iBoundary){
	if(mesh.boundaryMesh().types()[iBoundary] == "wall"){
		//Info << "Setting velocity at boundary patch back" << mesh.boundaryMesh().names()[iBoundary] << endl;
		U.boundaryField()[iBoundary] == U.boundaryField()[iBoundary] + Ugust.boundaryField()[iBoundary];
	}
}

Ugust == gustObj.getGustVelocity(time,U);
phiGust == (Ugust & mesh.Sf());

//Ensuring velocity is imposed on wall types patches
forAll(mesh.boundaryMesh(),iBoundary){
	if(mesh.boundaryMesh().types()[iBoundary] == "wall"){
		Info << "Applying velocity to boundary patch " << mesh.boundaryMesh().names()[iBoundary] << endl;
		U.boundaryField()[iBoundary] == U.boundaryField()[iBoundary] - Ugust.boundaryField()[iBoundary];
		phi.boundaryField()[iBoundary] = U.boundaryField()[iBoundary] & mesh.Sf().boundaryField()[iBoundary];
	}
}
