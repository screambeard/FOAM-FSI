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
    Foam::NPDSpectrumGustVelocityFvPatchVectorField

Description
    Foam::NPDSpectrumGustVelocityFvPatchVectorField

SourceFiles
    NPDSpectrumGustVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef NPDSpectrumGustVelocityFvPatchVectorField_H
#define NPDSpectrumGustVelocityFvPatchVectorField_H

#include "Random.H"
#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class NPDSpectrumGustFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class NPDSpectrumGustVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Reference value
        vectorField refValue_;

        //- Amplitude
        vectorField amplitude_;

        //- Direction
        vector direction_;

        //- Min frequency;
        scalar minFrequency_;

        //- Max frequency;
        scalar maxFrequency_;

        //- Number of segments
        label Nsegments_;

        //- Reference velocity (magnitude)
        scalar Uref_;

        //- random phase shifts
        List<scalar> phaseShifts_;

        //- Frequencies (size=Nsegments)
        List<scalar> frequencies_;

        //- Amplitudes (size=Nsegments)
        List<scalar> amplitudes_;

        //- Random generator
        Random rndGen_;

        //- Current time index
        label curTimeIndex_;

    // Private member functions

        //- Return current scale
        scalar currentAmplitude() const;

        //- Set frequences, amplitudes and phaseshifts
        void setSpectrum();

        scalar getAmplitude(scalar frequency) const;

public:

    //- Runtime type information
    TypeName("NPDSpectrumGustVelocity");


    // Constructors

        //- Construct from patch and internal field
        NPDSpectrumGustVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        NPDSpectrumGustVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given NPDSpectrumGustVelocityFvPatchVectorField
        //  onto a new patch
        NPDSpectrumGustVelocityFvPatchVectorField
        (
            const NPDSpectrumGustVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        NPDSpectrumGustVelocityFvPatchVectorField
        (
            const NPDSpectrumGustVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new NPDSpectrumGustVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        NPDSpectrumGustVelocityFvPatchVectorField
        (
            const NPDSpectrumGustVelocityFvPatchVectorField&,
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
                new NPDSpectrumGustVelocityFvPatchVectorField(*this, iF)
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

            //- Return amplitude
            vectorField amplitude() const
            {
                return amplitude_;
            }

            vectorField& amplitude()
            {
                return amplitude_;
            }

            //- Return frequency
            scalar maxFrequency() const
            {
                return maxFrequency_;
            }

            scalar& maxFrequency()
            {
                return maxFrequency_;
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
//#   include "NPDSpectrumGustVelocityFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
