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

\*---------------------------------------------------------------------------*/

#include "directionVaryingInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
#include "Tensor.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

vector directionVaryingInletOutletVelocityFvPatchVectorField::currentValue() const
{
    //Calculate angle of vector rotation on current time
    scalar time = dimensionedInternalField().mesh().time().value();
    scalar theta = amplitude_*Foam::sin(2*mathematicalConstant::pi*frequency_*time + phaseShift_);

    //Cross Product Matrix of direction vector
    Tensor<scalar> crPrMtrDir(vector::zero,vector::zero,vector::zero);
    crPrMtrDir.xy()=-direction_.z();
    crPrMtrDir.yx()=direction_.z();
    crPrMtrDir.xz()=direction_.y();
    crPrMtrDir.zx()=-direction_.y();
    crPrMtrDir.yz()=-direction_.x();
    crPrMtrDir.zy()=direction_.x();

    //Outer product of direction vectors
    Tensor<scalar> crPrDir = (direction_ * direction_);

    //Identy matrix with cos(theta) on diagonal
    Tensor<scalar> eye(vector::zero,vector::zero,vector::zero);
    eye.xx()=Foam::cos(theta);
    eye.yy()=Foam::cos(theta);
    eye.zz()=Foam::cos(theta);

    //Calculate new vector based on (http://en.wikipedia.org/wiki/Rotation_matrix#Basic_rotations)
    //This describes the rotation of a vector along an direction with a specified angle
    //rotation matrix = cos(theta)*I + sin(theta)*crPrMtrDir + (1- cos(theta))*crPrDir
    Tensor<scalar> rotM = eye + Foam::sin(theta)*crPrMtrDir + (1-Foam::cos(theta))*crPrDir;

    vector value = rotM & meanValue_;

    return value;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

directionVaryingInletOutletVelocityFvPatchVectorField::directionVaryingInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    meanValue_(vector::zero),
    direction_(vector(1,0,0)),
    frequency_(1),
    amplitude_(0),
    phaseShift_(0)
{
    this->refValue() = pTraits<vector>::zero;
    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 0.0;
}


directionVaryingInletOutletVelocityFvPatchVectorField::directionVaryingInletOutletVelocityFvPatchVectorField
(
    const directionVaryingInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    meanValue_(ptf.meanValue_),
    direction_(ptf.direction_),
    frequency_(ptf.frequency_),
    amplitude_(ptf.amplitude_),
    phaseShift_(ptf.phaseShift_)
{}


directionVaryingInletOutletVelocityFvPatchVectorField::directionVaryingInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    meanValue_(dict.lookup("meanValue")),
    direction_(dict.lookup("direction")),
    frequency_(readScalar(dict.lookup("frequency"))),
    amplitude_(readScalar(dict.lookup("amplitude"))),
    phaseShift_(readScalar(dict.lookup("phaseShift")))
{
    this->refValue() = vectorField(p.size(),currentValue());
    direction_ /= mag(direction_);
    amplitude_ *= mathematicalConstant::pi/180.0;
    phaseShift_ *= mathematicalConstant::pi/180.0;

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(this->refValue());
    }

    this->refGrad() = pTraits<vector>::zero;
    this->valueFraction() = 0.0;
}

directionVaryingInletOutletVelocityFvPatchVectorField::directionVaryingInletOutletVelocityFvPatchVectorField
(
    const directionVaryingInletOutletVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchField<vector>(ptf),
    phiName_(ptf.phiName_),
    meanValue_(ptf.meanValue_),
    direction_(ptf.direction_),
    frequency_(ptf.frequency_),
    amplitude_(ptf.amplitude_),
    phaseShift_(ptf.phaseShift_)
{}


directionVaryingInletOutletVelocityFvPatchVectorField::directionVaryingInletOutletVelocityFvPatchVectorField
(
    const directionVaryingInletOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    phiName_(ptf.phiName_),
    meanValue_(ptf.meanValue_),
    direction_(ptf.direction_),
    frequency_(ptf.frequency_),
    amplitude_(ptf.amplitude_),
    phaseShift_(ptf.phaseShift_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directionVaryingInletOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (!this->db().objectRegistry::found(phiName_))
    {
        // Flux not available, do not update
        mixedFvPatchField<vector>::updateCoeffs();

        return;
    }

    const scalarField& phip = this->lookupPatchField
    (
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    this->valueFraction() = 1.0 - pos(phip);

    this->refValue() = vectorField(phip.size(),currentValue());

    mixedFvPatchField<vector>::updateCoeffs();
}

void directionVaryingInletOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    if (phiName_ != "phi")
    {
        os.writeKeyword("phi")
            << phiName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("meanValue") << meanValue_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
    os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitude") << amplitude_*180.0/mathematicalConstant::pi << token::END_STATEMENT << nl;
    os.writeKeyword("phaseShift") << phaseShift_*180.0/mathematicalConstant::pi << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void directionVaryingInletOutletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& ptf
)
{
    fvPatchField<vector>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    directionVaryingInletOutletVelocityFvPatchVectorField
);


} // End namespace Foam

// ************************************************************************* //
