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
    rigidBody

Description
    Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    rigidBody.C

\*---------------------------------------------------------------------------*/

#ifndef rigidBody_H
#define rigidBody_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "simpleMatrix.H"
#include "septernion.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class rigidBody Declaration
\*---------------------------------------------------------------------------*/

class rigidBody
{
    // Private data
    const dictionary& dict_;
    const fvMesh& mesh_;
    const word name_;

    //Record time index to only return motion when timeIndex < mesh.timeIndex
    label timeIndex_;
    bool write_;

    //Patches of this body
    wordList patchNames_;
    labelList patchIDs_;

    //Solve properties
    word integrationScheme_;
    label interpolationRK4force_;

    //Body states
    scalarField bodyState_;
    scalarField bodyStatePrevIter_;
    scalarField bodyStateOld_;
    scalarField bodyStateOldOld_;
    scalarField bodyStateOldOldOld_;
    scalarField bodyState0_;

    //Body forces
    vector bodyForces_;
    vector bodyMoments_;
    vector bodyForcesOld_;
    vector bodyMomentsOld_;

    //Writing properties
    OFstream ofBody_;

    // Private Member Functions

        //- Check if interpolation option exists
        void checkInterpolationRK4force();

        //- Set patchIDs_ based on the patchNames
        void setPatchIDs();

        //- Do interpolation for forces or moments
        vector interpolateRK4term(const scalar rkCoeff,const vector valueOld, const vector value) const;

        //- Disallow default bitwise copy construct
        rigidBody(const rigidBody&);

    protected:

        const vector& calculateForces(const Field<vectorField>& patchForces);

        //- Interpolate function of RK forces
        vector interpolateRK4forces(const scalar rkCoeff) const;

        //- Interpolate function of RK moments
        vector interpolateRK4moments(const scalar rkCoeff) const;

        //- Return RK4 Coeffs
        SquareMatrix<scalar> getRK4Coeffs() const;

        //- Return BDF3 coeffs
        scalarField bdf3Coeffs() const;

        //- Return deltaT
        scalar deltaT_() const;

        //- Return deltaT0 or GREAT if timeIndex < 2 (first time step)
        scalar deltaT0_() const;

        scalarField& bodyState(){
            return bodyState_;}

        scalarField& bodyStateOld(){
            return bodyStateOld_;}

        scalarField& bodyStateOldOld(){
            return bodyStateOldOld_;}

        scalarField& bodyStateOldOldOld(){
            return bodyStateOldOldOld_;}

        scalarField& bodyState0(){
            return bodyState0_;}

        scalarField& bodyStatePrevIter(){
            return bodyStatePrevIter_;}

        vector& bodyForces(){
            return bodyForces_;}

        vector& bodyMoments(){
            return bodyMoments_;}

        OFstream& outputFile(){
            return ofBody_;}

        const dictionary& dict(){
            return dict_;}

        const fvMesh& mesh(){
            return mesh_;}

        const word& name(){
            return name_;}

        const word& integrationScheme(){
            return integrationScheme_;}

    public:

    //- Runtime type information
    TypeName("rigidBody");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            rigidBody,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh,
                const word name
            ),
            (dict, mesh, name)
        );


    // Selectors

        //- Return a pointer to the selected RigidMotion function
        static autoPtr<rigidBody> New
        (
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Constructors

        //- Construct from dict + mesh + name
        rigidBody(const dictionary& dict,const fvMesh& mesh,const word name);

        //- Create and return a clone
        //virtual autoPtr<rigidBody> clone() const = 0;

    // Destructor

        virtual ~rigidBody();

    // Member Functions

        //- Return patchIDs_;
        const labelList& getPatchIDs() const;

        //- Calculate forces for rigid body from patch forces
        virtual septernion calculateTransformation(const Field<vectorField>& patchForces) = 0;

        //- Update state variables for new time step
        virtual void update();

        //- Write out body properties if requested
        virtual void write() = 0;

        bool writeToFile() const;

        const label& localTimeIndex() const;
    };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
