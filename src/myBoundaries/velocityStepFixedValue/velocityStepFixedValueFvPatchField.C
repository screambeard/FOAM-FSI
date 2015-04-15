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

\*---------------------------------------------------------------------------*/

#include "velocityStepFixedValueFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar velocityStepFixedValueFvPatchField<Type>::currentScale() const
{
	scalar currentScale = 0.0;
	if(this->db().time().value()<(1.0/frequency_))
	{
		currentScale = 	1.0;
	}
    return currentScale;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
velocityStepFixedValueFvPatchField<Type>::velocityStepFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    refValue_(p.size()),
    amplitude_(p.size()),
    frequency_(0.0),
    curTimeIndex_(-1)
{}


template<class Type>
velocityStepFixedValueFvPatchField<Type>::velocityStepFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    refValue_("refValue", dict, p.size()),
    amplitude_("amplitude",dict, p.size()),
    frequency_(readScalar(dict.lookup("frequency"))),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(refValue_+amplitude_*currentScale());
    }
}


template<class Type>
velocityStepFixedValueFvPatchField<Type>::velocityStepFixedValueFvPatchField
(
    const velocityStepFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    amplitude_(ptf.amplitude_, mapper),
    frequency_(ptf.frequency_),
    curTimeIndex_(-1)
{}


template<class Type>
velocityStepFixedValueFvPatchField<Type>::velocityStepFixedValueFvPatchField
(
    const velocityStepFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    refValue_(ptf.refValue_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    curTimeIndex_(-1)
{}


template<class Type>
velocityStepFixedValueFvPatchField<Type>::velocityStepFixedValueFvPatchField
(
    const velocityStepFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    refValue_(ptf.refValue_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void velocityStepFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    refValue_.autoMap(m);
}


template<class Type>
void velocityStepFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const velocityStepFixedValueFvPatchField<Type>& tiptf =
        refCast<const velocityStepFixedValueFvPatchField<Type> >(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


template<class Type>
void velocityStepFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        patchField = refValue_+amplitude_*currentScale();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void velocityStepFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    amplitude_.writeEntry("amplitude", os);
//    os.writeKeyword("amplitude")
//        << amplitude_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
