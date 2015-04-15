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

Class
    Foam::directionVaryingInletOutletVelocityFvPatchVectorField

Description
    Foam::directionVaryingInletOutletVelocityFvPatchVectorField

SourceFiles
    directionVaryingInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef directionVaryingInletOutletVelocityFvPatchVectorField_H
#define directionVaryingInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class inletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class directionVaryingInletOutletVelocityFvPatchVectorField
:
    public mixedFvPatchVectorField
{

protected:

    // Protected data

        //- Name of flux field
        word phiName_;

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
    TypeName("directionVaryingInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        directionVaryingInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        directionVaryingInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given directionVaryingInletOutletVelocityFvPatchVectorField onto a new patch
        directionVaryingInletOutletVelocityFvPatchVectorField
        (
            const directionVaryingInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        directionVaryingInletOutletVelocityFvPatchVectorField
        (
            const directionVaryingInletOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new directionVaryingInletOutletVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        directionVaryingInletOutletVelocityFvPatchVectorField
        (
            const directionVaryingInletOutletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new directionVaryingInletOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
            }


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const fvPatchField<vector>& pvf);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
