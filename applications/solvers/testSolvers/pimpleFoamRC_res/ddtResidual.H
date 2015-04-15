volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) + (turbulence->divDevReff(U) & U) + fvc::grad(p);
scalarField magResU = mag( residual.internalField() );
scalar momentumResidual = std::sqrt( gSumSqr( magResU ) / mesh.globalData().nTotalCells());

// == Normalize residual == //
scalar Uav = ::sqrt( gSumSqr( mag(U.internalField()) )/ mesh.globalData().nTotalCells());
scalar momentumResidual1 =momentumResidual/(Uav/runTime.deltaT().value());
scalar dUdtAv = ::sqrt( gSumSqr( mag(fvc::ddt(U)().internalField() ) )/ mesh.globalData().nTotalCells() );
scalar momentumResidual2 = momentumResidual/dUdtAv;
scalar gradpAv = ::sqrt( gSumSqr( mag(fvc::grad(p)().internalField() ) )/ mesh.globalData().nTotalCells() );
scalar momentumResidual3 = momentumResidual/gradpAv;

Info << "dim|Udt|dUdt|gradP = " << momentumResidual << "|" << momentumResidual1 << "|" << momentumResidual2 << "|" << momentumResidual3 << endl;
Info << endl;
