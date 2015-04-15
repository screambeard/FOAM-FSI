const scalar& time = runTime.time().value();
const scalar& dt = runTime.deltaT().value();

Ugust == gustObj.getGustVelocity(time,U);

vectorField phiUgust = fvc::div(phi*Ugust)().internalField();


scalar freq = 5;
scalar yValue1 = 0.0;
scalar yValue0 = 0.0;
scalar c=1.0;//convective speed
if(time<1.0/freq){
	// Matlab function: u0(x0I)*(1-cos(2*pi*freq*(time-t0)))*0.5/(dx);
	scalar x1 = mesh.cellCentres()[rightCellIDs[0]].x();
	scalar x0 = mesh.cellCentres()[leftCellIDs[0]].x();
	yValue1 = 0.5*(1-Foam::cos(2*mathematicalConstant::pi*freq*time));
	yValue0 = 0.5*(1-Foam::cos(2*mathematicalConstant::pi*freq*(x0/c-x1/c-time)));
	//yValue = 1*time;
}

//Vortex characteristics
scalar C = 0.1;
scalar R = 0.1;
vector initPos(-0.1,0.0,0.5);//Check z-coordinate
vector curPos = initPos+vector(1,0,0)*time;//convective speed of 1 in x dir

//Setting left cells
gustSourceU == dimensionedVector("0",U.dimensions(),vector::zero);
forAll(rightCellIDs,iCell){
	label id = rightCellIDs[iCell];
	vector cc = mesh.cellCentres()[id];
	vector gustValue(0,0,0);

	//Vertical gust
	//gustValue.x() = 0;
	//gustValue.y() = yValue1;

	//Vortex gust
    scalar r=Foam::sqrt(pow(cc.x()-curPos.x(),2)+pow(cc.y()-curPos.y(),2));
	gustValue.x() = -(C*(cc.y()-curPos.y())/pow(R,2)) * Foam::exp(-0.5*pow(r/R,2));
	gustValue.y() = (C*(cc.x()-curPos.x())/pow(R,2)) * Foam::exp(-0.5*pow(r/R,2));

	gustSourceU.internalField()[id] = -gustValue;
}
volVectorField tmp1 = fvc::div(phi,gustSourceU);
forAll(leftCellIDs,iCell){
	label id = leftCellIDs[iCell];
	gustSource.internalField()[id] = -tmp1.internalField()[id];
}

//Setting right cells
gustSourceU == dimensionedVector("0",U.dimensions(),vector::zero);
forAll(leftCellIDs,iCell){
	label id = leftCellIDs[iCell];
	vector cc = mesh.cellCentres()[id];
	vector gustValue(0,0,0);

	//Vertical gust
	//gustValue.x() = 0;
	//gustValue.y() = yValue0;

	//Vortex gust
    scalar r=Foam::sqrt(pow(cc.x()-curPos.x(),2)+pow(cc.y()-curPos.y(),2));
	gustValue.x() = -(C*(cc.y()-curPos.y())/pow(R,2)) * Foam::exp(-0.5*pow(r/R,2));
	gustValue.y() = (C*(cc.x()-curPos.x())/pow(R,2)) * Foam::exp(-0.5*pow(r/R,2));


	gustSourceU.internalField()[id] = gustValue;
}
tmp1 == fvc::div(phi,gustSourceU);
forAll(rightCellIDs,iCell){
	label id = rightCellIDs[iCell];
	gustSource.internalField()[id] = -tmp1.internalField()[id];
}