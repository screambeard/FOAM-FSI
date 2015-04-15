    Info<< "Reading solid body properties\n" << endl;

    IOdictionary rigidBodyMotionDict
    (
        IOobject
        (
            "rigidBodyMotionDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    autoPtr<rigidBody> rbPtr
    (
        rigidBody::New
        (
            rigidBodyMotionDict,
            mesh
        )
    );

    rigidBody& rb = rbPtr();

    //- Save undisplaced points
    pointField undisplacedPoints_ = mesh.points();

    wordList farfieldPatches(rigidBodyMotionDict.lookup("farfieldPatches"));
    labelList farfieldPatchIDs(farfieldPatches.size(),-1);;
    forAll(farfieldPatches,ipatch)
    {
        label patchID(mesh.boundaryMesh().findPatchID(farfieldPatches[ipatch]));
        if(patchID < 0){
            FatalErrorIn("reading in farfieldPatches")
            << "patchName " << farfieldPatches[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            farfieldPatchIDs[ipatch] = patchID;
        }
    }
