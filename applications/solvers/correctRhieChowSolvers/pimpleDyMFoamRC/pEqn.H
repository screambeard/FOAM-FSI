HU = UEqn.H();
AU = UEqn.A();
U = HU/AU;

#include "calcPhi.H"

//adjustPhi(phi, U, p);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(1.0/fvc::interpolate(AU), p,"laplacian((1|A(U)),p)") == fvc::div(phi)
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
U -= (1.0/AU)*fvc::grad(p);
U.correctBoundaryConditions();
