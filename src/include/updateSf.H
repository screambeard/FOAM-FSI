{
    Sf.oldTime();

    const vectorField& allFaceAreas = mesh.faceAreas();

    Sf.internalField() =
        vectorField::subField(allFaceAreas, mesh.nInternalFaces());

    const fvPatchList& patches = mesh.boundary();

    forAll (patches, patchI)
    {
        Sf.boundaryField()[patchI] = 
            patches[patchI].patchSlice(allFaceAreas);
    }
}

