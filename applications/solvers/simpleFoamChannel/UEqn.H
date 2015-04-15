    // Solve the Momentum equation

	dimensionedVector gradP
    (
        "gradP",
        dimensionSet(0, 1, -2, 0, 0),
        vector(1.0,0,0)
    );

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
      ==
      gradP
    );

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
