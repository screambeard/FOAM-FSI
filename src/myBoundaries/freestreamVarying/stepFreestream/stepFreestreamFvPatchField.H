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
    Foam::stepFreestreamFvPatchField

Description
    Foam::stepFreestreamFvPatchField

SourceFiles
    stepFreestreamFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef stepFreestreamFvPatchField_H
#define stepFreestreamFvPatchField_H

#include "freestreamFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class freestreamFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class stepFreestreamFvPatchField
:
    public freestreamFvPatchField<Type>
{

    // Private data

        //- Amplitude
        Field<Type> endValue_;

        //- time to reach endValue
        scalar dT_;

        //- Current time index
        label curTimeIndex_;

        //- Time when constructed: increase will hapen from time0_ to time0_+dT_
        scalar time0_;

    // Private member functions

        //- Return current scale
        scalar currentScale() const;


public:

    //- Runtime type information
    TypeName("stepFreestream");


    // Constructors

        //- Construct from patch and internal field
        stepFreestreamFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        stepFreestreamFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given stepFreestreamFvPatchField onto a new patch
        stepFreestreamFvPatchField
        (
            const stepFreestreamFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        stepFreestreamFvPatchField
        (
            const stepFreestreamFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new stepFreestreamFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        stepFreestreamFvPatchField
        (
            const stepFreestreamFvPatchField<Type>&,
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
                new stepFreestreamFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Return defining fields

            const Field<Type>& freestreamValue() const
            {
                return this->refValue();
            }

            Field<Type>& freestreamValue()
            {
                return this->refValue();
            }

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
#   include "stepFreestreamFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
