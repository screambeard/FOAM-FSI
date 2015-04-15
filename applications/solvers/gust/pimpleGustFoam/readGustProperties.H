Info<< "Reading gustProperties\n" << endl;
IOdictionary gustProperties
(
    IOobject
    (
        "gustProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

word gustType(gustProperties.lookup("gustType"));

// === gust object === //
autoPtr<gust> gustPtr
(
    gust::New
    (
        gustType,
        gustProperties,
        mesh
    )
);
gust& gustObj = gustPtr();

// Create gust flux field
surfaceScalarField phiGust
(
    IOobject
    (
        "phiGust",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    U.dimensions() * mesh.Sf().dimensions(),
    fixedValueFvPatchVectorField::typeName
);
phiGust == dimensionedScalar("phiGust",U.dimensions() * mesh.Sf().dimensions(),0.0);

// Create dummyField used to determine ddt Coefficient for the DGCL
surfaceScalarField phiDummy1
(
    IOobject
    (
        "phiDummy1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mesh.phi().dimensions(),
    fixedValueFvPatchVectorField::typeName
);
phiDummy1 == dimensionedScalar("phiDummy1",phiDummy1.dimensions(),1.0);

// Create dummyField used to determine ddt Coefficient for the DGCL
surfaceScalarField phiDummy0
(
    IOobject
    (
        "phiDummy0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    mesh.phi().dimensions(),
    fixedValueFvPatchVectorField::typeName
);
phiDummy0 == dimensionedScalar("phiDummy0",phiDummy0.dimensions(),0.0);

//Field to hold artificial change in volume due to phiGust (DGCL)
volScalarField VolumeGust
(
    IOobject
    (
        "VolumeGust",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimVolume,
    zeroGradientFvPatchScalarField::typeName
);
VolumeGust.internalField() = mesh.V();
VolumeGust.oldTime().internalField() = mesh.V();
VolumeGust.oldTime().oldTime().internalField() = mesh.V();

//Used in ddt scheme since change in Volume is time depentent. Holds ratio between VolumeGust and mesh.V for all used times
volScalarField volumeRatio
(
    IOobject
    (
        "volumeRatio",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimless,
    zeroGradientFvPatchScalarField::typeName
);
volumeRatio.internalField() = 1.0;
volumeRatio.oldTime().internalField() = 1.0;
volumeRatio.oldTime().oldTime().internalField() = 1.0;
