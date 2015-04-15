volScalarField rAU = 1.0/UEqn.A();
U = rAU*UEqn.H();
phi = (fvc::interpolate(U) & mesh.Sf());
// ddtPhiCorr does not work.  HJ, 20/Nov/2013

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAU, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    eqnResidual = pEqn.solve().initialResidual();
    //Only use first of nonOrth corrector initial residuals. This residual is used for the PISO loops, not the non-orth
    if(nonOrth==0){
        maxResidual = max(eqnResidual, maxResidual);
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}

//Only relax if not last corrector
if(eqnResidual > innerConvergence && corr < nCorr){
    //Info << "Doing Relexation on p" << endl;
	p.relax();
}

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);
U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
