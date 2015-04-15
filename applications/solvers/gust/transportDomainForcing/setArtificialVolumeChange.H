// ======== Set correct volume change (mesh.phi()) ======== //
mesh.phi();
mesh.setPhi() == phiDummy1;
scalar phi1_0 = fvc::meshPhi(U)()[0];
mesh.setPhi() == phiDummy0;
surfaceScalarField phi0  = fvc::meshPhi(U);
scalar phi0_0 = phi0[0];
scalar sweptVolumeCoeff = phi1_0 - phi0_0;//ddtScheme dependent

//Calculate the required swept volume at t(n+1) for the DGCL based on phiGust
mesh.setPhi() == (-phiGust - phi0)/sweptVolumeCoeff;

//Calculate addedVolume based on mesh.phi() and calculate volumeRatio
volScalarField addedVolume = fvc::div(mesh.phi())*runTime.deltaT();
addedVolume.internalField() *= mesh.V();//fvc operator devides by mesh.V()
VolumeGust.internalField() = VolumeGust.internalField() + addedVolume.internalField();
volumeRatio.internalField() = VolumeGust.internalField()/mesh.V();