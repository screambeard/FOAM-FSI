    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    U.oldTime();

    Info<< "Creating surfaceField Ugust\n" << endl;
    surfaceVectorField Ugust
    (
        IOobject
        (
            "UgustFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        U.dimensions(),
        fixedValueFvPatchVectorField::typeName
    );
    Ugust == dimensionedVector("UgustFace",U.dimensions(),vector::zero);
    Ugust.oldTime() == dimensionedVector("UgustFace",U.dimensions(),vector::zero);

    volVectorField Utotal
    (
        IOobject
        (
            "Utotal",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U+fvc::average(Ugust)
    );
    Utotal.write();

#   include "createPhi.H"
    phi.oldTime();

    //Creating field to store old face areas needed voor ddtPhiCorr
#   include "createSf.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    //UEqn H and A holders
    Info<< "Reading field AU if present\n" << endl;
    volScalarField AU
    (
        IOobject
        (
            "AU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        1.0/runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    );
    
    Info<< "Reading field HU if present\n" << endl;
    volVectorField HU
    (
        IOobject
        (
            "HU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        U.dimensions()/runTime.deltaT().dimensions(),
        zeroGradientFvPatchVectorField::typeName
    );

    volScalarField Udummy
    (
        IOobject
        (
            "Udummy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0,0,0,0,0,0,0),
        fixedValueFvPatchVectorField::typeName
    );
    Udummy == 1.0;
    Udummy.oldTime() == 1.0;

    //Turbulence model
    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );