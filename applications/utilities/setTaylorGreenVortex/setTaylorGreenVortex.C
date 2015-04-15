/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    applyBoundaryLayer

Description
    Apply a simplified boundary-layer model to the velocity and
    turbulence fields based on the 1/7th power-law.

    The uniform boundary-layer thickness is either provided via the -ybl option
    or calculated as the average of the distance to the wall scaled with
    the thickness coefficient supplied via the option -Cbl.  If both options
    are provided -ybl is used.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDist.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "initContinuityErrs.H"

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );
    U.oldTime();

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    #include "createPhi.H"
    phi.oldTime();

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

    // =======================================================================//
    // Calculate U (u,v) based on:  u = sin(x).*cos(y);
    //                              v = -cos(x).*sin(y);
    //                              p = 0.25*(cos(2*x)+cos(2*y));
    // Assumed that 0 < x,y < 2*pi
    // =======================================================================//


    const vectorField& cc(mesh.cellCentres());
    forAll(U, celli)
    {
        U[celli].x() = Foam::sin(cc[celli].x())*Foam::cos(cc[celli].y());
        U[celli].y() = -Foam::cos(cc[celli].x())*Foam::sin(cc[celli].y());
        p[celli] = 0.25*(Foam::cos(2*cc[celli].x()) + Foam::cos(2*cc[celli].y()));
    }
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();

    //Calculate phi
    phi = fvc::interpolate(U)&mesh.Sf();

    //Creating field to store old face areas needed voor ddtPhiCorr
#   include "createSf.H"

    bool solvePressure=false;
    label nOuter = 10;
    for(int i=0;i<nOuter;i++){

        fvVectorMatrix UEqn
        (
          fvm::div(phi, U)
          - fvm::laplacian(nu, U)
          //+ turbulence->divDevReff(U)
        );

        volVectorField HU = UEqn.H();
        volScalarField AU = UEqn.A();
        U = HU/AU;

        //Calculate non-divergence free flux
        phi = fvc::interpolate(HU)/fvc::interpolate(AU) & mesh.Sf();

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

        //Build pEqn
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1.0/fvc::interpolate(AU), p,"laplacian((1|A(U)),p)") == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        if(solvePressure){
            pEqn.solve();
        }
        phi -= pEqn.flux();

        //if(i!=nOuter-1){
            //Reset U again
            forAll(U, celli)
            {
                U[celli].x() = Foam::sin(cc[celli].x())*Foam::cos(cc[celli].y());
                U[celli].y() = -Foam::cos(cc[celli].x())*Foam::sin(cc[celli].y());
            }
        //}else{
        //    U -= fvc::grad(p)/AU;
        //}
        U.correctBoundaryConditions();
    }

    //Write out fields
    phi.write();
    U.write();
    p.write();

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
