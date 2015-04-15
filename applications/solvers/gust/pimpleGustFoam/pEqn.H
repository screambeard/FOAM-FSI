HU = UEqn().H();
AU = UEqn().A();
U = HU/AU;

if (nCorr <= 1)
{
    UEqn.clear();
}

// ==== Two important steps compard to icoFoamRC: make phi absolute and setDdtPhiCorr using old Sf (manually for EulerRC) ==== //
//Make absolute phi before using in ddtPhiCorr
phi = phi - phiGust;

#include "setDdtPhiCorr.H"
phi = ((fvc::interpolate(HU)/fvc::interpolate(AU)) & mesh.Sf());
phi += ddtPhiCorrLocal;

forAll(phi.boundaryField(),iPatch)
{
	if(!phi.boundaryField()[iPatch].coupled())
	{
		phi.boundaryField()[iPatch] = U.boundaryField()[iPatch] & mesh.Sf().boundaryField()[iPatch];
	}
}

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
        outCorr == nOuterCorr - 1
     && corr == nCorr - 1
     && nonOrth == nNonOrthCorr
    )
    {
        pEqn.solve(mesh.solutionDict().solver("pFinal"));
    }
    else
    {
        pEqn.solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector except for last corrector
if (outCorr != nOuterCorr-1)
{
    p.relax();
}

//Make relative phi
phi = phi + phiGust;

U -= (1.0/AU)*fvc::grad(p);
U.correctBoundaryConditions();
