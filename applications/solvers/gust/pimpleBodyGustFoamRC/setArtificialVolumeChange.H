// ======== Set correct volume change (mesh.phi()) ======== //
//Save mesh movement fluxes and faceSweptVolumes
const surfaceScalarField phiMesh(fvc::meshPhi(U));//mesh flux
const surfaceScalarField meshphi(mesh.phi());
mesh.setPhi().oldTime() == mesh.setPhi().oldTime() + meshPhiGust.oldTime();
const surfaceScalarField phiMeshPlus(fvc::meshPhi(U));
mesh.setPhi().oldTime() == mesh.setPhi().oldTime() - meshPhiGust.oldTime();

//Get coefficient
mesh.setPhi() == phiDummy1;
scalar phi1_0 = fvc::meshPhi(U)()[0];
mesh.setPhi() == phiDummy0;
scalar phi0_0 = fvc::meshPhi(U)()[0];
scalar sweptVolumeCoeff = phi1_0 - phi0_0;//ddtScheme dependent

mesh.setPhi() == meshphi;//set it back to mesh phi

//Calculate the required swept volume at t(n+1) for the DGCL based on phiGust
//const surfaceScalarField meshphiGust((phiGust - phi0)/sweptVolumeCoeff);
meshPhiGust == (-phiGust + phiMesh - phiMeshPlus)/sweptVolumeCoeff;

//Calculate addedVolume based on mesh.phi() and calculate volumeRatio
volScalarField addedVolume = fvc::div(meshPhiGust)*runTime.deltaT();
addedVolume.internalField() *= mesh.V();//fvc operator devides by mesh.V()

VolumeGust.internalField() = VolumeGust.internalField() + addedVolume.internalField();
volumeRatio.internalField() = (mesh.V() + VolumeGust.internalField())/mesh.V();