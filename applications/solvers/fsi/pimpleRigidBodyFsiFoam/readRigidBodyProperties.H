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
