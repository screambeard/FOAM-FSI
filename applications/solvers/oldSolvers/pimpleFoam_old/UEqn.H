// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

if(!lastIter){
    UEqn().relax();
}

solve(UEqn() == -fvc::grad(p));
