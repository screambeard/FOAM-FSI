{
    // PIMPLE loop
	//Outer correction loops needed for correct order convergence -> pimplefoam like
	label oCorr=1;
	scalar relativeResidual = 1;
	bool converged=false;
    bool lastIter=false;
	do
	{
		Info << "outer iteration: " << oCorr << endl;
        maxResidual = 0;//reset maxResidual for new outerloop

		fvc::makeRelative(phi, U);
		p.storePrevIter();
		U.storePrevIter();

		// phi RELATIVE
		tmp<fvVectorMatrix> UEqn
		(
        	fvm::ddt(U)
      		+ fvm::div(phi, U)
      		+ turbulence->divDevReff(U)
		);

		if(relaxOnLastIter || !lastIter){
    		UEqn().relax();
		}

		solve(UEqn() == -fvc::grad(p));

		UEqn() =
		(
		    fvm::ddt(U)
		  + fvm::div(phi, U)
		  + turbulence->divDevReff(U)
		);

        // --- PISO loop
        label corr=0;
		do
		{
			HU = UEqn().H();
		    AU = UEqn().A();
	        U = HU/AU;

	    #	include "calcPhi.H"

			//In adjustPhi already request p.needReference(). If not nothing is done. If yes than first do makeRelative, do stuff and than makeAbsolute
            //adjustPhi(phi, U, p);

		    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
		    {
		        fvScalarMatrix pEqn
		        (
		            fvm::laplacian(1.0/fvc::interpolate(AU), p,"laplacian((1|A(U)),p)") == fvc::div(phi)
		        );

		        pEqn.setReference(pRefCell, pRefValue);
		        eqnResidual = pEqn.solve().initialResidual();
    			maxResidual = max(eqnResidual, maxResidual);
				corr++;//Update corrector number

		        if (nonOrth == nNonOrthCorr)
		        {
		            phi -= pEqn.flux();
		        }
		    }

			//Only relax if not last corrector
			if(relaxOnLastIter || (eqnResidual > innerConvergence && corr != nCorr) || !lastIter){
            	p.relax();
			}

	//#       include "movingMeshContinuityErrs.H"
            fvc::makeRelative(phi, U);

		    U -= 1.0/AU*fvc::grad(p);
		    U.correctBoundaryConditions();
		}while(eqnResidual > innerConvergence && corr < nCorr);



		turbulence->correct();
		fvc::makeAbsolute(phi, U);

		scalar residualPressure = gSumMag( p.internalField() - p.prevIter().internalField() ) / gSumMag( p.internalField() );
    	scalar residualVelocity = gSumMag( U.internalField() - U.prevIter().internalField() ) / gSumMag( U.internalField() );
    	relativeResidual = max( residualPressure, residualVelocity );
    	Info << "residualPressure = " << residualPressure << ", residualVelocity = " << residualVelocity << endl;

		//Check if this was the last iter and set converged to true
		if(lastIter){
			converged = true;
		}
		//Check convergence and set lastIter to true
		if(relativeResidual < outerConvergence || oCorr == nOuterCorr){
			lastIter = true;
		}

		oCorr++;
	}while(!converged);

}
