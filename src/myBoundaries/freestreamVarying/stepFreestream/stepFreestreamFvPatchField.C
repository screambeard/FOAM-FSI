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

#include "stepFreestreamFvPatchField.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar stepFreestreamFvPatchField<Type>::currentScale() const
{
    scalar currentScale = 0.0;
    if(this->db().time().value()<time0_)
    {
        currentScale = 0.0;
    }
    else if(this->db().time().value()<(time0_+dT_))
    {
        currentScale =  0.5*(1.0 - cos(2*mathematicalConstant::pi*(1/dT_/2)*(this->db().time().value()-time0_) ));
    }
    else
    {
        currentScale = 1.0;
    }
    return currentScale;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
stepFreestreamFvPatchField<Type>::stepFreestreamFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    freestreamFvPatchField<Type>(p, iF),
    endValue_(p.size()),
    dT_(0.0),
    curTimeIndex_(-1),
    time0_(this->db().time().value())
{}


template<class Type>
stepFreestreamFvPatchField<Type>::stepFreestreamFvPatchField
(
    const stepFreestreamFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    freestreamFvPatchField<Type>(ptf, p, iF, mapper),
    endValue_(ptf.endValue_, mapper),
    dT_(ptf.dT_),
    curTimeIndex_(-1),
    time0_(ptf.time0_)
{}


template<class Type>
stepFreestreamFvPatchField<Type>::stepFreestreamFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    freestreamFvPatchField<Type>(p, iF),
    endValue_("endValue",dict, p.size()),
    dT_(readScalar(dict.lookup("dT"))),
    curTimeIndex_(-1),
    time0_(readScalar(dict.lookup("startTime")))
{
    freestreamValue() = Field<Type>("freestreamValue", dict, p.size());

    if (dict.found("value"))
    {
        freestreamFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        freestreamFvPatchField<Type>::operator==(freestreamValue() + (endValue_ - freestreamValue())*currentScale());
    }

    if (dict.found("phi"))
    {
        dict.lookup("phi") >> this->phiName_;
    }
}


template<class Type>
stepFreestreamFvPatchField<Type>::stepFreestreamFvPatchField
(
    const stepFreestreamFvPatchField<Type>& ptf
)
:
    freestreamFvPatchField<Type>(ptf),
    endValue_(ptf.endValue_),
    dT_(ptf.dT_),
    curTimeIndex_(-1),
    time0_(ptf.time0_)
{}


template<class Type>
stepFreestreamFvPatchField<Type>::stepFreestreamFvPatchField
(
    const stepFreestreamFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    freestreamFvPatchField<Type>(ptf, iF),
    endValue_(ptf.endValue_),
    dT_(ptf.dT_),
    curTimeIndex_(-1),
    time0_(ptf.time0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void stepFreestreamFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        this->refValue() = freestreamValue() + (endValue_ - freestreamValue())*currentScale();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    freestreamFvPatchField<Type>::updateCoeffs();
}

template<class Type>
void stepFreestreamFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (this->phiName_ != "phi")
    {
        os.writeKeyword("phi")
            << this->phiName_ << token::END_STATEMENT << nl;
    }
    freestreamValue().writeEntry("freestreamValue", os);
    endValue_.writeEntry("endValue", os);
    os.writeKeyword("dT") << dT_ << token::END_STATEMENT << nl;
    os.writeKeyword("startTime") << time0_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
