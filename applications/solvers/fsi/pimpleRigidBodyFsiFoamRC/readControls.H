#   include "readTimeControls.H"
#   include "readPIMPLEControls.H"
#   include "readPIMPLELoopControls.H"

    bool correctPhi = false;
    if (pimple.found("correctPhi"))
    {
        correctPhi = Switch(pimple.lookup("correctPhi"));
    }

    bool checkMeshCourantNo = false;
    if (pimple.found("checkMeshCourantNo"))
    {
        checkMeshCourantNo = Switch(pimple.lookup("checkMeshCourantNo"));
    }

    bool relaxOnLastIter = pimple.lookupOrDefault("relaxOnLastIteration",false);
    Info << "relaxOnLastIteration = " << relaxOnLastIter << endl;
