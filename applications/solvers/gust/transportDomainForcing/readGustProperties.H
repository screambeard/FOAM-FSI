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

volVectorField gustSource
(
    IOobject
    (
        "gustSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    phi.dimensions()*U.dimensions()/dimVolume,
    zeroGradientFvPatchVectorField::typeName
);
gustSource==dimensionedVector("gustSource",phi.dimensions()*U.dimensions()/dimVolume,vector::zero);

volVectorField gustSourceU
(
    IOobject
    (
        "gustSourceU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    U.dimensions(),
    zeroGradientFvPatchVectorField::typeName
);
gustSourceU==dimensionedVector("gustSourceU",U.dimensions(),vector::zero);

// Base forcing cells on set
word forceSetName="forceSet";
faceSet forceFaceSet(mesh,"forceFaceSet",IOobject::MUST_READ);

volScalarField gamma
(
    IOobject
    (
        "gamma",
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0),
    zeroGradientFvPatchScalarField::typeName
);

//Determine left and right cells
vector planeNormal(1,0,0);
labelList leftCellIDs(forceFaceSet.size(),-1);
labelList rightCellIDs(forceFaceSet.size(),-1);
const labelList& nei = mesh.neighbour();
const labelList& own = mesh.owner();

label counter = 0;
forAll(forceFaceSet,iFace)
{
    label faceID = forceFaceSet.toc()[iFace];
    vector ownCC = mesh.cellCentres()[own[faceID]];
    vector neiCC = mesh.cellCentres()[nei[faceID]];
    vector faceC = mesh.faceCentres()[faceID];
    scalar ownInProduct = ((ownCC - faceC) & planeNormal);
    scalar neiInProduct = ((ownCC - faceC) & planeNormal);
    if(ownInProduct>0){
        rightCellIDs[counter] = own[faceID];
        leftCellIDs[counter] = nei[faceID];

        gamma.internalField()[own[faceID]] = 1;
        gamma.internalField()[nei[faceID]] = -1;
    }else{
        rightCellIDs[counter] = nei[faceID];
        leftCellIDs[counter] = own[faceID];

        gamma.internalField()[own[faceID]] = -1;
        gamma.internalField()[nei[faceID]] = 1;
    }
    counter++;
}