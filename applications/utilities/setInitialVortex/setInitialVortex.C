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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("freestreamVelocity");
    argList::validArgs.append("position");
    argList::validArgs.append("C");
    argList::validArgs.append("R");


#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "initContinuityErrs.H"

    vector freestreamVelocity(IStringStream(args.additionalArgs()[0])());
    vector vortexPosition(IStringStream(args.additionalArgs()[1])());
    scalar vortexStrength(readScalar(IStringStream(args.additionalArgs()[2])()));
    scalar vortexRadius(readScalar(IStringStream(args.additionalArgs()[3])()));

    Info << "freestreamVelocity = " << freestreamVelocity << endl;
    Info << "vortexPosition = " << vortexPosition << endl;
    Info << "vortexStrength = " << vortexStrength << endl;
    Info << "vortexRadius = " << vortexRadius << endl;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

#   include "createPhi.H"

    surfaceVectorField UFace
    (
        IOobject
        (
            "UFace",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(U)
    );

#   include "continuityErrs.H"
    // =======================================================================//
    //Calculate U (u,v) based on:   u = Uinf - C(y-yc)/R^2 * exp(-0.5*(r/R)^2)
    //                              v = C(x -xc)/R^2 * exp(-0.5*(r/R)^2)
    // =======================================================================//

    bool cellEvaluation = true;
    if(cellEvaluation){
        Info << "Evaluating vortex at cellcentres and calculating phi by interpolation" << endl;
        const vectorField& cc(mesh.cellCentres());
        scalarField radius = mag(cc - vortexPosition);
        forAll(U, celli)
        {
            U[celli].x() = U[celli].x() - vortexStrength*(cc[celli].y() - vortexPosition.y())/sqr(vortexRadius) * Foam::exp(-0.5*sqr(radius[celli]/vortexRadius));
            U[celli].y() = U[celli].y() + vortexStrength*(cc[celli].x() - vortexPosition.x())/sqr(vortexRadius) * Foam::exp(-0.5*sqr(radius[celli]/vortexRadius));
            p[celli] = p[celli] - 0.5*sqr(vortexStrength/vortexRadius) * Foam::exp(-sqr(radius[celli]/vortexRadius));
        }

        phi = fvc::interpolate(U) & mesh.Sf();
    }else{
        Info << "Evaluating vortex at faces and reconstructing from there" << endl;
        const vectorField& fc(mesh.faceCentres());
        scalarField radius = mag(fc - vortexPosition);

        //InternalField
        forAll(UFace.internalField(),facei){
            UFace.internalField()[facei].x() = freestreamVelocity.x() - vortexStrength*(fc[facei].y() - vortexPosition.y())/sqr(vortexRadius) * Foam::exp(-0.5*sqr(radius[facei]/vortexRadius));
            UFace.internalField()[facei].y() = vortexStrength*(fc[facei].x() - vortexPosition.x())/sqr(vortexRadius) * Foam::exp(-0.5*sqr(radius[facei]/vortexRadius));
        }

        //BoundaryField
        forAll(UFace.boundaryField(),iBoundary){
            forAll(UFace.boundaryField()[iBoundary],facei){
                label faceID = mesh.boundaryMesh()[iBoundary].start() + facei;
                UFace.boundaryField()[iBoundary][facei].x() = freestreamVelocity.x() - vortexStrength*(fc[faceID].y() - vortexPosition.y())/sqr(vortexRadius) * Foam::exp(-0.5*sqr(radius[faceID]/vortexRadius));
                UFace.boundaryField()[iBoundary][facei].y() = vortexStrength*(fc[faceID].x() - vortexPosition.x())/sqr(vortexRadius) * Foam::exp(-0.5*sqr(radius[faceID]/vortexRadius));
            }
        }
        phi = (UFace & mesh.Sf());
        U = fvc::reconstruct(phi);
    }

    volScalarField divPhi
    (
        IOobject
        (
            "divPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(phi)
    );
    divPhi.write();

    U.write();
    phi.write();
    p.write();


#   include "continuityErrs.H"

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
