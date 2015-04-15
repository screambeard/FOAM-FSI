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

Class
    twoDoFRigidBody2Masses

Description
    Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    twoDoFRigidBody2Masses.C

\*---------------------------------------------------------------------------*/

#ifndef twoDoFRigidBody2Masses_H
#define twoDoFRigidBody2Masses_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "rigidBody.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"
#include "Time.H"
#include "simpleMatrix.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class twoDoFRigidBody2Masses Declaration
\*---------------------------------------------------------------------------*/

class twoDoFRigidBody2Masses
:
    public rigidBody
{
    //Body properties
    scalar k_;
    scalar c_;
    scalar ktheta_;
    scalar ctheta_;
    scalar mtheta_;
    scalar m_;
    scalar Irc_;
    scalar Lcg2rc_;
    vector rc_;
    scalar thetaGeom_;//In radians

    // Private Member Functions

        //- Writing function
        void writeState(const scalar& time);

        //- Write header
        void writeHeader();

        //- Return only force in y direction
        scalar getForce(const vector forces) const;

        //- Return only force in y direction
        scalar getMoment(const vector moments) const;

        //- Calculate moment around z
        void calculateMoments(const Field<vectorField>& patchForces);

        //- Function for damping coefficient
        scalar calculateThetaDamping(const scalar& theta) const;

        //- Calculate motion
        scalar calculatePosition();

        //- Read initial state from dict if present (state = [x,y,theta,dxdt,dydt,dthetadt])
        void readInitialStateIfPresent(const dictionary& dict);

        //- Disallow default bitwise copy construct
        twoDoFRigidBody2Masses(const twoDoFRigidBody2Masses&);

    public:

    //- Runtime type information
    TypeName("twoDoFRigidBody2Masses");

    //Public data

    // Constructors

        //- Construct from dict + mesh + name
        twoDoFRigidBody2Masses(const dictionary& dict,const fvMesh& mesh,const word name);

    // Destructor

        virtual ~twoDoFRigidBody2Masses();

    // Member Functions

        //- Calculate forces for rigid body from patch forces
        virtual septernion calculateTransformation(const Field<vectorField>& patchForces);

        //- Write out body properties if requested
        virtual void write();
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
