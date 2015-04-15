// Reading/calculating Sf field for ddtPhiCorr
IOobject VolumesHeader
(
    "Volumes",
    runTime.timeName(),
    mesh,
    IOobject::READ_IF_PRESENT,
    IOobject::AUTO_WRITE
);

volScalarField Volumes
(
    VolumesHeader,
    mesh,
    dimensionedScalar("0", dimVolume, 0.0)
);

if(!VolumesHeader.headerOk())
{
    Volumes.internalField() = mesh.V();
}
