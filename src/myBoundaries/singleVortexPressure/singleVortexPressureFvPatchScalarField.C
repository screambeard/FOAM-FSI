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

#include "singleVortexPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalarField Foam::singleVortexPressureFvPatchScalarField::currentAmplitude() const
{
    scalarField vortexPressure(patch_.size(),0.0);

    //Calculate radius and all face values (internal & boundary)
    scalarField radius = mag(patch_.Cf() - position_);
    vortexPressure = - 0.5*sqr(strength_/R_) * Foam::exp(-sqr(radius/R_));

    return vortexPressure;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singleVortexPressureFvPatchScalarField::singleVortexPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    patch_(p),
    UName_("U"),
    p0_(p.size()),
    convSpeed_(vector::zero),
    R_(0),
    strength_(0),
    position_(vector::zero),
    curTimeIndex_(-1)
{}


Foam::singleVortexPressureFvPatchScalarField::singleVortexPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    patch_(p),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    p0_("p0", dict, p.size()),
    convSpeed_(dict.lookup("convectionSpeed")),
    R_(readScalar(dict.lookup("radius"))),
    strength_(readScalar(dict.lookup("strength"))),
    position_(dict.lookup("position")),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);//Needs something extra
    }
}


Foam::singleVortexPressureFvPatchScalarField::singleVortexPressureFvPatchScalarField
(
    const singleVortexPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    patch_(p),
    UName_(ptf.UName_),
    p0_(ptf.p0_, mapper),
    convSpeed_(ptf.convSpeed_),
    R_(ptf.R_),
    strength_(ptf.strength_),
    position_(ptf.position_),
    curTimeIndex_(-1)
{}


Foam::singleVortexPressureFvPatchScalarField::singleVortexPressureFvPatchScalarField
(
    const singleVortexPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    patch_(ptf.patch_),
    UName_(ptf.UName_),
    p0_(ptf.p0_),
    convSpeed_(ptf.convSpeed_),
    R_(ptf.R_),
    strength_(ptf.strength_),
    position_(ptf.position_),
    curTimeIndex_(-1)
{}


Foam::singleVortexPressureFvPatchScalarField::singleVortexPressureFvPatchScalarField
(
    const singleVortexPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    patch_(ptf.patch_),
    UName_(ptf.UName_),
    p0_(ptf.p0_),
    convSpeed_(ptf.convSpeed_),
    R_(ptf.R_),
    strength_(ptf.strength_),
    position_(ptf.position_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singleVortexPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}


void Foam::singleVortexPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const singleVortexPressureFvPatchScalarField& tiptf =
        refCast<const singleVortexPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::singleVortexPressureFvPatchScalarField::updateCoeffs(const vectorField& Up)
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        scalarField& patchField = *this;

        position_ += convSpeed_*db().time().deltaT().value();
        Info << "Vortex position = " << position_ << endl;

        patchField = p0_+currentAmplitude();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::singleVortexPressureFvPatchScalarField::updateCoeffs()
{
    updateCoeffs(patch().lookupPatchField<volVectorField, vector>(UName_));
}


void Foam::singleVortexPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    p0_.writeEntry("p0", os);
    os.writeKeyword("convectionSpeed")
        << convSpeed_ << token::END_STATEMENT << nl;
    os.writeKeyword("radius")
        << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("strength")
        << strength_ << token::END_STATEMENT << nl;
    os.writeKeyword("position")
        << position_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        singleVortexPressureFvPatchScalarField
    );
}

// ************************************************************************* //
