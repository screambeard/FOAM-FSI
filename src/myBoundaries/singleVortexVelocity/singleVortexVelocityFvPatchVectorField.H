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
    Foam::singleVortexVelocityFvPatchVectorField

Description
    Foam::singleVortexVelocityFvPatchVectorField

SourceFiles
    singleVortexVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef singleVortexVelocityFvPatchVectorField_H
#define singleVortexVelocityFvPatchVectorField_H

#include "Random.H"
#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NPDSpectrumGustFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class singleVortexVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Patch reference
        const fvPatch& patch_;

        //- Reference value
        vectorField refValue_;

        //- ConvectionSpeed
        vector convSpeed_;

        //- Core radius
        scalar R_;

        //- Circulation of vortex
        scalar strength_;

        //- Initial position
        vector position_;

        //- Current time index
        label curTimeIndex_;

    // Private member functions

        //- Return current scale
        vectorField currentAmplitude() const;

public:

    //- Runtime type information
    TypeName("singleVortexVelocity");


    // Constructors

        //- Construct from patch and internal field
        singleVortexVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        singleVortexVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given singleVortexVelocityFvPatchVectorField
        //  onto a new patch
        singleVortexVelocityFvPatchVectorField
        (
            const singleVortexVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        singleVortexVelocityFvPatchVectorField
        (
            const singleVortexVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new singleVortexVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        singleVortexVelocityFvPatchVectorField
        (
            const singleVortexVelocityFvPatchVectorField&,
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
                new singleVortexVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const vectorField& refValue() const
            {
                return refValue_;
            }

            //- Return reference to the ref value to allow adjustment
            vectorField& refValue()
            {
                return refValue_;
            }

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
//#   include "singleVortexVelocityFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
