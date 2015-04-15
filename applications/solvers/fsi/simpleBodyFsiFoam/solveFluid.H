{
    //SIMPLE loop
	tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

	solve(UEqn() == -fvc::grad(p));

	//Pressure equation
	p.boundaryField().updateCoeffs();

	AU = UEqn().A();
    U = UEqn().H()/AU;
    UEqn.clear();
    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, p);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1.0/fvc::interpolate(AU), p,"laplacian((1|A(U)),p)") == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

#	include "continuityErrs.H"

	// Explicitly relax pressure for momentum corrector
	p.relax();

    U -= 1.0/AU*fvc::grad(p);
    U.correctBoundaryConditions();
}
turbulence->correct();