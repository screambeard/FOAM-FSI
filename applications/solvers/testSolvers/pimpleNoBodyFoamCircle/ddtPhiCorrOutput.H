IOobject ddtIOobject
    (
        "ddtPhiCorr(" + AU.name() + ',' + U.name() + ',' + phi.name() + ')',
        mesh.time().timeName(),
        mesh
    );

scalar deltaT = mesh.time().deltaT().value();

scalar deltaT0 = mesh.time().deltaT0().value();
if (U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex() || U.oldTime().oldTime().timeIndex() < 0)
{
    deltaT0 = GREAT;
}

scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
scalar coefft0  = coefft + coefft00;

dimensionedScalar rDeltaT = 1.0/mesh.time().deltaT();

//Set coefficient for excluding certain boundaries from ddtPhiCorr
surfaceScalarField ddtPhiCoeff
(
    IOobject
    (
        "ddtPhiCoeff",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensioned<scalar>("1", dimless, 1.0)
);

//Set coefficients at boundary which have fixedValues to 0
forAll(U.boundaryField(), patchI)
{
    if
    (
        U.boundaryField()[patchI].fixesValue()
     //|| isA<symmetryFvPatchVectorField>(U.boundaryField()[patchI])
     //|| isA<slipFvPatchVectorField>(U.boundaryField()[patchI])
    )
    {
        ddtPhiCoeff.boundaryField()[patchI] = 0.0;
    }
}

//Set old velocity interpolated to new mesh
surfaceVectorField U0f = fvc::interpolate(U.oldTime());

//Set old old velocity interpolated to new mesh
surfaceVectorField U00f = fvc::interpolate(U.oldTime().oldTime());

//Get Sf with old times recored
const surfaceVectorField& Sf =
    mesh.objectRegistry::lookupObject<surfaceVectorField>
    (
        "Sf"
    );

//Corr test6
surfaceScalarField phi0 = (phi.oldTime() - (U0f&mesh.Sf()));
surfaceScalarField phi00 = (phi.oldTime().oldTime() - (U00f&mesh.Sf()));

surfaceScalarField ddtPhiCorr6
(
    ddtIOobject,
    rDeltaT*ddtPhiCoeff
   *(
        coefft0*phi0
      - coefft00*phi00
    )/fvc::interpolate(AU)
);

//Corr test7
phi0 = (phi.oldTime() - (U0f&Sf.oldTime()));
phi00 = (phi.oldTime().oldTime() - (U00f&Sf.oldTime().oldTime()));

surfaceScalarField ddtPhiCorr7
(
    ddtIOobject,
    rDeltaT*ddtPhiCoeff
   *(
        coefft0*phi0
      - coefft00*phi00
    )/fvc::interpolate(AU)
);

//Corr test8
phi0 = (phi.oldTime() - phi.oldTime());
phi00 = (phi.oldTime().oldTime() - phi.oldTime().oldTime());

surfaceScalarField ddtPhiCorr8
(
    ddtIOobject,
    rDeltaT*ddtPhiCoeff
   *(
        coefft0*phi0
      - coefft00*phi00
    )/fvc::interpolate(AU)
);

scalar norm6=sqrt(sum(sqr(ddtPhiCorr6))).value();
scalar norm7=sqrt(sum(sqr(ddtPhiCorr7))).value();
scalar norm8=sqrt(sum(sqr(ddtPhiCorr8))).value();

diffPhiCorrOutput << runTime.time().value() << "\t" << norm7 << "\t" << norm6 << "\t" << norm8 << endl;
