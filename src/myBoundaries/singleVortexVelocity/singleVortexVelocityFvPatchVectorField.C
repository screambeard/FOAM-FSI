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

#include "singleVortexVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

vectorField singleVortexVelocityFvPatchVectorField::currentAmplitude() const
{
	vectorField vortexVelocities(patch_.size(),vector::zero);
    /*const scalar& pi = mathematicalConstant::pi;
    vectorField relativePosition = patch_.Cf() - position_;
    scalarField radius = mag(relativePosition);
    scalarField angle = atan2(relativePosition.component(1),relativePosition.component(0));
    scalarField radialVelocity = circulation_/2/pi*(radius/(sqr(radiusVortexCore_) + sqr(radius)));

    forAll(vortexVelocities,iFace){
        vortexVelocities[iFace].x() = -radialVelocity[iFace]*Foam::sin(angle[iFace]);
        vortexVelocities[iFace].y() = radialVelocity[iFace]*Foam::cos(angle[iFace]);
    }*/

    //Calculate radius and all face values (internal & boundary)
    scalarField radius = mag(patch_.Cf() - position_);
    vortexVelocities.replace(0,-strength_*(patch_.Cf().component(1) - position_.y())/sqr(R_) * Foam::exp(-0.5*sqr(radius/R_)));
    vortexVelocities.replace(1,strength_*(patch_.Cf().component(0) - position_.x())/sqr(R_) * Foam::exp(-0.5*sqr(radius/R_)));

    return vortexVelocities;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

singleVortexVelocityFvPatchVectorField::singleVortexVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    patch_(p),
    refValue_(p.size()),
    convSpeed_(vector::zero),
    R_(0),
    strength_(0),
    position_(vector::zero),
    curTimeIndex_(-1)
{}


singleVortexVelocityFvPatchVectorField::singleVortexVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    patch_(p),
    refValue_("refValue", dict, p.size()),
    convSpeed_(dict.lookup("convectionSpeed")),
    R_(readScalar(dict.lookup("radius"))),
    strength_(readScalar(dict.lookup("strength"))),
    position_(dict.lookup("position")),
    curTimeIndex_(-1)
{
    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==(vectorField("value", dict, p.size()));
    }
    else
    {
        fixedValueFvPatchVectorField::operator==(refValue_+currentAmplitude());
    }
}

singleVortexVelocityFvPatchVectorField::singleVortexVelocityFvPatchVectorField
(
    const singleVortexVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    patch_(p),
    refValue_(ptf.refValue_, mapper),
    convSpeed_(ptf.convSpeed_),
    R_(ptf.R_),
    strength_(ptf.strength_),
    position_(ptf.position_),
    curTimeIndex_(-1)
{}

singleVortexVelocityFvPatchVectorField::singleVortexVelocityFvPatchVectorField
(
    const singleVortexVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    patch_(ptf.patch_),
    refValue_(ptf.refValue_),
    convSpeed_(ptf.convSpeed_),
    R_(ptf.R_),
    strength_(ptf.strength_),
    position_(ptf.position_),
    curTimeIndex_(-1)
{}

singleVortexVelocityFvPatchVectorField::singleVortexVelocityFvPatchVectorField
(
    const singleVortexVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    patch_(ptf.patch_),
    refValue_(ptf.refValue_),
    convSpeed_(ptf.convSpeed_),
    R_(ptf.R_),
    strength_(ptf.strength_),
    position_(ptf.position_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void singleVortexVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    refValue_.autoMap(m);
}

void singleVortexVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const singleVortexVelocityFvPatchVectorField& tiptf =
        refCast<const singleVortexVelocityFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}

void singleVortexVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        vectorField& patchField = *this;

        position_ += convSpeed_*db().time().deltaT().value();
        Info << "Vortex position = " << position_ << endl;

        patchField = refValue_+currentAmplitude();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void singleVortexVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    refValue_.writeEntry("refValue", os);
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

makePatchTypeField
(
    fvPatchVectorField,
    singleVortexVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
