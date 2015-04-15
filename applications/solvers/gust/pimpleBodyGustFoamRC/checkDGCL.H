{
	//Calculate the change in artificial volume 
	const volScalarField& Vr = volumeRatio;
	//volScalarField conserve = -(fvc::div(fvc::meshPhi(U)));
	volScalarField divPhiMeshGust = -(fvc::div(fvc::meshPhi(U)-phiGust));
	volScalarField conserve = divPhiMeshGust;
	volScalarField ddtV = fvc::ddt(Vr, Udummy);
	conserve.internalField() += ddtV;

	//Check if div(phiGust) and change in volume is the same
	label nNonConservedCells = 0;
	forAll(conserve.internalField(),iCell){
		if(conserve.internalField()[iCell] > 1e-10/runTime.deltaT().value()){
			nNonConservedCells++;
			const point& cc = mesh.cellCentres()[iCell];
			Info << iCell << " " << cc << " :ddt(V) + div(phiMeshGust) = " << ddtV.internalField()[iCell] << " + " << divPhiMeshGust.internalField()[iCell] << " = " << conserve.internalField()[iCell] << " > " << 1e-12/runTime.deltaT().value() << endl;
		}
	}
	if(nNonConservedCells > 0){
		Info << "DGCL NOT CONSERVED: " << nNonConservedCells << " cells are not conserving" << endl;
	}
}