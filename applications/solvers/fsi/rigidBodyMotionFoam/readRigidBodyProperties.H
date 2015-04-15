    Info<< "Reading solid body properties\n" << endl;

    //Create rigid body
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

    // Read in forcing options
    bool forcing = rigidBodyMotionDict.lookupOrDefault("forcing",false);

    vector prescribedForces = vector::zero;
    scalar prescribedFrequency = 0;
    if(forcing){
        Info << "Add forcing to system" << endl;
        prescribedForces = vector(rigidBodyMotionDict.lookup("prescribedForces"));
        prescribedFrequency = readScalar(rigidBodyMotionDict.lookup("prescribedFrequency"));
    }
