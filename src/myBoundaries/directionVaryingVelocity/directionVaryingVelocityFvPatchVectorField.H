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
    Foam::directionVaryingVelocityFvPatchVectorField

Description
    Foam::directionVaryingVelocityFvPatchVectorField

SourceFiles
    directionVaryingVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef directionVaryingVelocityFvPatchVectorField_H
#define directionVaryingVelocityFvPatchVectorField_H

#include "Random.H"
#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NPDSpectrumGustFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class directionVaryingVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- mean velocity value
        vector meanValue_;

        //- oscillating direction
        vector direction_;

        //- frequency
        scalar frequency_;

        //- amplitude in degrees
        scalar amplitude_;

        //- phase-shift in degrees
        scalar phaseShift_;

        //- Return function for current value
        vector currentValue() const;

public:

    //- Runtime type information
    TypeName("directionVaryingVelocity");


    // Constructors

        //- Construct from patch and internal field
        directionVaryingVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        directionVaryingVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given directionVaryingVelocityFvPatchVectorField
        //  onto a new patch
        directionVaryingVelocityFvPatchVectorField
        (
            const directionVaryingVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        directionVaryingVelocityFvPatchVectorField
        (
            const directionVaryingVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new directionVaryingVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        directionVaryingVelocityFvPatchVectorField
        (
            const directionVaryingVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new directionVaryingVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "directionVaryingVelocityFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
