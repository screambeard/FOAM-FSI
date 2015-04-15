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
    Foam::travelingFixedValueFvPatchField

Description
    Foam::travelingFixedValueFvPatchField

SourceFiles
    travelingFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef travelingFixedValueFvPatchField_H
#define travelingFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class travelingFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class travelingFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        const fvPatch& patch_;

        //- Reference value
        Field<Type> refValue_;

        //- Amplitude
        Field<Type> amplitude_;

        //- Frequency
        scalar frequency_;

        //- ConvectionSpeed
        vector convSpeed_;

        //- Gust length
        scalar gustLength_;

        //- Current time index
        label curTimeIndex_;

        //- Minium patch x
        scalar minPatchX_;

    // Private member functions

        //- Return current field
        Field<Type> currentField() const;

        //- Function setting minimum x coordinate of patch
        void setMinPatchX();


public:

    //- Runtime type information
    TypeName("travelingFixedValue");


    // Constructors

        //- Construct from patch and internal field
        travelingFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        travelingFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given travelingFixedValueFvPatchField
        //  onto a new patch
        travelingFixedValueFvPatchField
        (
            const travelingFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        travelingFixedValueFvPatchField
        (
            const travelingFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new travelingFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        travelingFixedValueFvPatchField
        (
            const travelingFixedValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new travelingFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the ref value
            const Field<Type>& refValue() const
            {
                return refValue_;
            }

            //- Return reference to the ref value to allow adjustment
            Field<Type>& refValue()
            {
                return refValue_;
            }

            //- Return amplitude
            scalar amplitude() const
            {
                return amplitude_;
            }

            scalar& amplitude()
            {
                return amplitude_;
            }

            //- Return frequency
            scalar frequency() const
            {
                return frequency_;
            }

            scalar& frequency()
            {
                return frequency_;
            }

            //- Return frequency
            scalar convectionSpeed() const
            {
                return convSpeed_;
            }

            scalar& convectionSpeed()
            {
                return convSpeed_;
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
                const fvPatchField<Type>&,
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

#ifdef NoRepository
#   include "travelingFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
