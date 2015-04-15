//U.correctBoundaryConditions();
forAll(mesh.boundaryMesh(),ipatch){
    bool correctB = true;
    forAll(farfieldPatchIDs,fpatch){
        if(ipatch == farfieldPatchIDs[fpatch]){
            correctB = false;
        }
    }

    if(correctB){
        U.boundaryField()[ipatch].updateCoeffs();
    }
}
