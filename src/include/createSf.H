// Reading/calculating Sf field for ddtPhiCorr
IOobject SfHeader
(
    "Sf",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

surfaceVectorField Sf
(
    SfHeader,
    mesh,
    dimensionedVector("0", dimArea, vector::zero)
);

if(!SfHeader.headerOk())
{
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

