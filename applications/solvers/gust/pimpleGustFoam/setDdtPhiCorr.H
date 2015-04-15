dimensionedScalar rDeltaT = 1.0/mesh.time().deltaT();

IOobject ddtIOobject
(
    "ddtPhiCorr(" + AU.name() + ',' + U.name() + ',' + phi.name() + ')',
    mesh.time().timeName(),
    mesh
);

// === Set boundaries correct of U === //
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

forAll (U.boundaryField(), patchI)
{
    ddtPhiCoeff.boundaryField()[patchI] = 0.0;
}

surfaceScalarField ddtPhiCorrLocal
(
	ddtIOobject,
	mesh,
	dimensioned<scalar>(word("1"),dimensionSet(0,3,-1,0,0,0,0),1.0)
);

if(ddtScheme==word("EulerRC"))
{
	tmp<surfaceScalarField> phiCorr =
		phi.oldTime()*fvc::interpolate(volumeRatio.oldTime()) - (fvc::interpolate(U.oldTime()*volumeRatio.oldTime()) & mesh.Sf());//TEMP: Everything should be old including Sf (Paper tukovic about Rhie-Chow interpolation)
	//surfaceScalarField phiCorr =
	//    phi.oldTime() - (fvc::interpolate(U.oldTime()) & mesh.Sf());//TEMP: Everything should be old including Sf (Paper tukovic about Rhie-Chow interpolation)

	// === make sure fvc::interpolate is performed on 1.0/rA seperately === //
//	tmp<surfaceScalarField> ddtPhiCorrLocal = tmp<surfaceScalarField>

	ddtPhiCorrLocal = tmp<surfaceScalarField>
	(
		new surfaceScalarField
		(
		    ddtIOobject,
		    rDeltaT*ddtPhiCoeff
		   *1.0/fvc::interpolate(AU)*phiCorr
		)
	);
}
else if(ddtScheme==word("backwardRC"))
{
	// === Set coefficients (copied from backwardRC) === //
	scalar deltaT = mesh.time().deltaT().value();
	scalar deltaT0 = 0.0;
	if(U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex() || U.oldTime().oldTime().timeIndex() < 0)
	{
		deltaT0 = GREAT;
	}
	else
	{
		deltaT0 = mesh.time().deltaT0().value();
	}

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

	// === make sure fvc::interpolate is performed on 1.0/rA seperately === //
	//tmp<surfaceScalarField> ddtPhiCorrLocal = tmp<surfaceScalarField>
	ddtPhiCorrLocal = tmp<surfaceScalarField>
	(
		new surfaceScalarField
		(
            ddtIOobject,
            rDeltaT*ddtPhiCoeff
           *(
                (
		           	(
		               coefft0*fvc::interpolate(volumeRatio.oldTime())*phi.oldTime()
		             - coefft00*fvc::interpolate(volumeRatio.oldTime().oldTime())*phi.oldTime().oldTime()
		            )
		          - (       
                        ((coefft0*fvc::interpolate(U.oldTime()*volumeRatio.oldTime())) & mesh.Sf())
                      - ((coefft00*fvc::interpolate(U.oldTime().oldTime()*volumeRatio.oldTime().oldTime())) & mesh.Sf())
		            )
				)/fvc::interpolate(AU)
            )
        )
    );
}
