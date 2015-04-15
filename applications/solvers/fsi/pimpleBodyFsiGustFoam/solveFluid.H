{
    // PIMPLE loop
	//Outer correction loops needed for correct order convergence -> pimplefoam like
	label oCorr=0;
	scalar relativeResidual = 1;
	do
	{
		Info << "outer iteration: " << oCorr+1 << endl;
        maxResidual = 0;//reset maxResidual for new outerloop

		// Make the fluxes relative
        fvc::makeRelative(phi, U);
        phi = phi + phiGust;

		p.storePrevIter();
		U.storePrevIter();

		fvVectorMatrix UEqn
		(
        	fvm::ddt(volumeRatio,U)
      		+ fvm::div(phi, U)
      		+ turbulence->divDevReff(U)
		);

		UEqn.relax();

		solve(UEqn == -fvc::grad(p));

		// --- PISO loop
        label corr=0;
		do
		{
		    HU = UEqn.H();
		    AU = UEqn.A();
	        U = HU/AU;

        #	include "calcPhiGust.H"

		#	include "adjustPhiGust.H"

		    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
		    {
		        fvScalarMatrix pEqn
		        (
		            fvm::laplacian(1.0/fvc::interpolate(AU), p,"laplacian((1|A(U)),p)") == fvc::div(phi)
		        );

		        pEqn.setReference(pRefCell, pRefValue);
		        eqnResidual = pEqn.solve().initialResidual();//Solve equation
    			maxResidual = max(eqnResidual, maxResidual);//Update residual
				corr++;//Update corrector number

		        if (nonOrth == nNonOrthCorr)
		        {
		            phi -= pEqn.flux();
		        }
		    }

			//Only relax if not last corrector
			if(eqnResidual > innerConvergence && corr < nCorr){
            	p.relax();
			}

	#       include "movingMeshContinuityErrs.H"

		    // Make the fluxes relative
		    fvc::makeRelative(phi, U);
		    phi = phi + phiGust;

		    U -= 1.0/AU*fvc::grad(p);
		    U.correctBoundaryConditions();


		}while(eqnResidual > innerConvergence && corr < nCorr);

		turbulence->correct();

		// Make the fluxes absolute
        phi = phi - phiGust;//gust
        fvc::makeAbsolute(phi, U);

		scalar residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / gSumMag( p.internalField() );
    	scalar residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / gSumMag( U.internalField() );
    	relativeResidual = max( residualPressure, residualVelocity );
    	Info << "residualPressure = " << residualPressure << endl;
    	Info << "residualVelocity = " << residualVelocity << endl;
		
		oCorr++;

	}while((relativeResidual > outerConvergence || oCorr < 2) && oCorr < nOuterCorr);

}
