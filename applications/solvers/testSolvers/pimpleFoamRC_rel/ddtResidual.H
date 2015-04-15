volVectorField residual = fvc::ddt( U ) + fvc::div( phi, U ) + (turbulence->divDevReff(U) & U) + fvc::grad(p);
scalarField magResU = mag( residual.internalField() );
scalar momentumResidual = ::sqrt( sum( sqr( magResU ) ) / mesh.nCells() );

// == Normalize residual == //
//scalar Uav = ::sqrt( sum( sqr( mag(U.internalField()) ) )/ mesh.nCells() );
//momentumResidual /= Uav/runTime.deltaT().value();
//scalar dUdtAv = ::sqrt( sum( sqr( mag(fvc::ddt(U)().internalField() ) ) )/ mesh.nCells() );
//momentumResidual /= dUdtAv;
scalar gradpAv = ::sqrt( sum( sqr( mag(fvc::grad(p)().internalField() ) ) )/ mesh.nCells() );
momentumResidual /= gradpAv;

bool convergence = momentumResidual <= outerConvergence;

labelList convergenceList( Pstream::nProcs(), 0 );
convergenceList[Pstream::myProcNo()] = convergence;
reduce( convergenceList, sumOp<labelList>() );

int minIter = 2;
convergence = min( convergenceList ) && oCorr >= minIter - 1;

Info << "root mean square residual norm = " << momentumResidual << ", tolerance = " << outerConvergence;
Info << ", iteration = " << oCorr + 1;
Info << ", convergence = ";

if ( convergence )
  Info << "true";
else
  Info << "false";

Info << endl;
