HU = UEqn().H();
AU = UEqn().A();
U = HU/AU;

if (nCorr <= 1)
{
    UEqn.clear();
}

#include "calcPhiGust.H"

#include "adjustPhiGust.H"

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(1.0/fvc::interpolate(AU), p,"laplacian((1|A(U)),p)") == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    
    if
    (
        ocorr == nOuterCorr-1
     && corr == nCorr-1
     && nonOrth == nNonOrthCorr
    )
    {
        pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
    }
    else
    {
        pEqn.solve(mesh.solutionDict().solver(p.name()));
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector except for last corrector
if (ocorr != nOuterCorr-1)
{
    p.relax();
}

fvc::makeRelative(phi,U);
phi = phi + phiGust;

U -= (1.0/AU)*fvc::grad(p);
U.correctBoundaryConditions();