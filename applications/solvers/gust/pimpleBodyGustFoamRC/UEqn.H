tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(volumeRatio,U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

if (ocorr != nOuterCorr-1)
{
    UEqn().relax();
}

if (momentumPredictor)
{
    solve(UEqn() == -fvc::grad(p));
}