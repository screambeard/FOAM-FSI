/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    wallShearStress

Description
    Calculates and writes the wall shear stress, for the specified times.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.append("freestreamVelocity");

    timeSelector::addOptions();
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    vector Uinf(IStringStream(args.additionalArgs()[0])());

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        #include "createFields.H"

        volSymmTensorField Reff(RASModel->devReff());

        volVectorField wallShearStress
        (
            IOobject
            (
                "wallShearStress",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "wallShearStress",
                Reff.dimensions(),
                vector::zero
            )
        );

        volScalarField skinFriction
        (
            IOobject
            (
                "Cf",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(0)
        );

        forAll(wallShearStress.boundaryField(), patchi)
        {
            wallShearStress.boundaryField()[patchi] =
            (
                mesh.Sf().boundaryField()[patchi]
                /mesh.magSf().boundaryField()[patchi]
            ) & Reff.boundaryField()[patchi];
        }


        //Determine correct tangential direction for all faces
        surfaceVectorField nf = mesh.Sf()/mesh.magSf();
        surfaceVectorField tf = nf;
        forAll(nf.boundaryField(),patchi){

            tf.boundaryField()[patchi].replace(0,nf.boundaryField()[patchi].component(1));
            tf.boundaryField()[patchi].replace(1,-nf.boundaryField()[patchi].component(0));
            
            scalarField tfSign = sign(tf.boundaryField()[patchi] & Uinf);
            tf.boundaryField()[patchi] *= tfSign;

            scalarField tauSign = sign(tf.boundaryField()[patchi] & wallShearStress.boundaryField()[patchi]);

            skinFriction.boundaryField()[patchi] = mag(wallShearStress.boundaryField()[patchi])*tauSign;
            //skinFriction.boundaryField()[patchi] = (wallShearStress.boundaryField()[patchi] & tf.boundaryField()[patchi]);
        }
        skinFriction /= (0.5*sqr(mag(Uinf)));

        skinFriction.write();

    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
