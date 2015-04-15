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

#include "travelingFixedValueFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Field<Type> travelingFixedValueFvPatchField<Type>::currentField() const
{
    Field<Type> fieldCur(patch_.size());
    scalar time = this->db().time().value();

    //TEMP: assume convSpeed is in x direction only
    scalar xGustFront = convSpeed_.x()*time+minPatchX_;
    Info << "Gust front at " << xGustFront << endl;

    forAll(patch_.Cf(),iFace){
        const vector& position(patch_.Cf()[iFace]);
        scalar currentNormX = position.x()-xGustFront+gustLength_;

        if(currentNormX >= 0 && currentNormX <= gustLength_){
            fieldCur[iFace] = amplitude_[iFace]*0.5*(1-Foam::cos(2*mathematicalConstant::pi/gustLength_*currentNormX));
        }else{
            fieldCur[iFace] = amplitude_[iFace]*0.0;
        }
    }

    return fieldCur;
}

template<class Type>
void travelingFixedValueFvPatchField<Type>::setMinPatchX(){
    const pointField& points = patch_.boundaryMesh().mesh().points();
    const labelList& mp = patch_.patch().meshPoints();
    const pointField patchPoints(points,mp);
    minPatchX_ = Foam::min(patchPoints.component(0));
    Info << "min x of patch = " << minPatchX_ << endl;
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
travelingFixedValueFvPatchField<Type>::travelingFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    patch_(p),
    refValue_(p.size()),
    amplitude_(p.size()),
    frequency_(0.0),
    convSpeed_(vector::zero),
    gustLength_(mag(convSpeed_)/frequency_),
    curTimeIndex_(-1),
    minPatchX_(0)
{
    setMinPatchX();
}


template<class Type>
travelingFixedValueFvPatchField<Type>::travelingFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    patch_(p),
    refValue_("refValue", dict, p.size()),
    amplitude_("amplitude",dict, p.size()),
    frequency_(readScalar(dict.lookup("frequency"))),
    convSpeed_(dict.lookup("convectionSpeed")),
    gustLength_(mag(convSpeed_)/frequency_),
    curTimeIndex_(-1),
    minPatchX_(0)
{
    setMinPatchX();

    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(refValue_+currentField());
    }
}


template<class Type>
travelingFixedValueFvPatchField<Type>::travelingFixedValueFvPatchField
(
    const travelingFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    patch_(p),
    refValue_(ptf.refValue_, mapper),
    amplitude_(ptf.amplitude_, mapper),
    frequency_(ptf.frequency_),
    convSpeed_(ptf.convSpeed_),
    gustLength_(mag(convSpeed_)/frequency_),
    curTimeIndex_(-1)
{
    setMinPatchX();
}


template<class Type>
travelingFixedValueFvPatchField<Type>::travelingFixedValueFvPatchField
(
    const travelingFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    patch_(ptf.patch_),
    refValue_(ptf.refValue_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    convSpeed_(ptf.convSpeed_),
    gustLength_(mag(convSpeed_)/frequency_),
    curTimeIndex_(-1)
{
    setMinPatchX();
}


template<class Type>
travelingFixedValueFvPatchField<Type>::travelingFixedValueFvPatchField
(
    const travelingFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    patch_(ptf.patch_),
    refValue_(ptf.refValue_),
    amplitude_(ptf.amplitude_),
    frequency_(ptf.frequency_),
    convSpeed_(ptf.convSpeed_),
    gustLength_(mag(convSpeed_)/frequency_),
    curTimeIndex_(-1)
{
    setMinPatchX();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void travelingFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    refValue_.autoMap(m);
}


template<class Type>
void travelingFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const travelingFixedValueFvPatchField<Type>& tiptf =
        refCast<const travelingFixedValueFvPatchField<Type> >(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


template<class Type>
void travelingFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        patchField = refValue_+currentField();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void travelingFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    amplitude_.writeEntry("amplitude", os);
    os.writeKeyword("frequency")
        << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("convectionSpeed")
        << convSpeed_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
