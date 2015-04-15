    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    if(!lastIter){
        //Info << "Doing Relexation on U" << endl;
        UEqn.relax();
    }
    solve(UEqn == -fvc::grad(p));
