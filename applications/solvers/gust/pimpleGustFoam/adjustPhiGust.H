if (p.needReference())
{
    p.boundaryField().updateCoeffs();
    scalar massIn = 0.0;
    scalar fixedMassOut = 0.0;
    scalar adjustableMassOut = 0.0;

    //make relative to phiGust do this after p.boundaryField().updateCoeffs()!
    phi = phi + phiGust;

    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  HJ, 13/Feb/2009
    if (phi.mesh().moving())
    {
        fvc::makeRelative(phi, U);
    }
    forAll (phi.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        const fvsPatchScalarField& phip = phi.boundaryField()[patchi];

        // Bug fix: All coupled patches should also be unaffected
        // HJ, 12/Feb/2010
        if (!phip.coupled())
        {
            if
            (
                Up.fixesValue()
             && !isA<inletOutletFvPatchVectorField>(Up)
            )
            {
                forAll (phip, i)
                {
                    if (phip[i] < 0.0)
                    {
                        massIn -= phip[i];
                    }
                    else
                    {
                        fixedMassOut += phip[i];
                    }
                }
            }
            else
            {
                forAll(phip, i)
                {
                    if (phip[i] < 0.0)
                    {
                        massIn -= phip[i];
                    }
                    else
                    {
                        adjustableMassOut += phip[i];
                    }
                }
            }
        }
    }

    //HR 16.03.10: Bug fix for moving meshes. Changing domain volume needs
    // to be taken into account.
    if (phi.mesh().moving())
    {
        dimensionedScalar Vdiff =
            sum(phi.mesh().V()*volumeRatio.internalField()) - sum(phi.mesh().V0()*volumeRatio.oldTime().internalField());

        fixedMassOut += Vdiff.value()/phi.time().deltaT().value();
    }
    else
    {
        dimensionedScalar Vdiff =
            sum(phi.mesh().V()*volumeRatio.internalField()) - sum(phi.mesh().V()*volumeRatio.oldTime().internalField());

        fixedMassOut += Vdiff.value()/phi.time().deltaT().value();
    }

    reduce(massIn, sumOp<scalar>());
    reduce(fixedMassOut, sumOp<scalar>());
    reduce(adjustableMassOut, sumOp<scalar>());

    scalar massCorr = 1.0;
    static const scalar closedDomainTol =
        debug::tolerances("closedDomainTol", 1e-10);

    if (mag(adjustableMassOut) > SMALL)
    {
        massCorr = (massIn - fixedMassOut)/adjustableMassOut;
    }
    else if
    (
        mag(fixedMassOut - massIn)
      > closedDomainTol*Foam::max(1.0, mag(massIn))
    )
    {
        // Cannot adjust
        FatalErrorIn
        (
            "adjustPhi\n"
            "(\n"
            "    surfaceScalarField& phi,\n"
            "    const volVectorField& U,\n"
            "    const volScalarField& p\n"
            ")"
        )   << "Continuity error cannot be removed by adjusting the"
            " outflow.\nPlease check the velocity boundary conditions"
            " and/or run potentialFoam to initialise the outflow." << nl
            << "Specified mass inflow   : " << massIn << nl
            << "Specified mass outflow  : " << fixedMassOut << nl
            << "Difference              : "
            << mag(fixedMassOut - massIn) << nl
            << "Adjustable mass outflow : " << adjustableMassOut << nl
            << exit(FatalError);
    }

    if (fvMesh::debug)
    {
        Info<< "bool Foam::adjustPhi(...) massIn: " << massIn
            << " fixedMassOut: " << fixedMassOut
            << " adjustableMassOut: " << adjustableMassOut
            << " mass corr: " << massCorr
            << endl;
    }
    
    forAll (phi.boundaryField(), patchi)
    {
        const fvPatchVectorField& Up = U.boundaryField()[patchi];
        fvsPatchScalarField& phip = phi.boundaryField()[patchi];

        // Bug fix: All coupled patches should also be unaffected
        // HJ, 12/Feb/2010
        if (!phip.coupled())
        {
            if
            (
                !Up.fixesValue()
             || isA<inletOutletFvPatchVectorField>(Up)
            )
            {
                forAll (phip, i)
                {
                    if (phip[i] > 0.0)
                    {
                        phip[i] *= massCorr;
                    }
                }
            }
        }
    }
    
    // If the mesh is moving, adjustment needs to be calculated on
    // relative fluxes.  Now reverting to absolute fluxes.  HJ, 13/Feb/2009
    if (phi.mesh().moving())
    {
        fvc::makeAbsolute(phi, U);
    }
    
    //makeAbsolute after adjustphi
    phi = phi - phiGust;
}
// ************************************************************************* //
