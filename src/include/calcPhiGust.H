{
    phi = (fvc::interpolate(HU)/fvc::interpolate(AU)) & mesh.Sf();
    
    forAll(phi.boundaryField(), patchI)
    {
        if (!phi.boundaryField()[patchI].coupled())
        {
            phi.boundaryField()[patchI] =
            (
                U.boundaryField()[patchI] 
              & mesh.Sf().boundaryField()[patchI]
            );
        }
    }

#   include "setDdtPhiCorrGust.H"
    phi += ddtPhiCorrLocal;
}

