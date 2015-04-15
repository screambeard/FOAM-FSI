tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(volumeRatio,U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

//Relax UEqn
if (outCorr == nOuterCorr-1)
{
    UEqn().relax(1);
}
else
{
    UEqn().relax();
}

//Solve equation with final solver if last iteration
if (outCorr == nOuterCorr-1)
{
    solve(UEqn() == -fvc::grad(p), mesh.solutionDict().solver("UFinal"));
}
else
{
    solve(UEqn() == -fvc::grad(p));
}
