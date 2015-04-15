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

#include "directionVaryingVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"
#include "Tensor.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
vector directionVaryingVelocityFvPatchVectorField::currentValue() const
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

directionVaryingVelocityFvPatchVectorField::directionVaryingVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    meanValue_(vector::zero),
	direction_(vector(1,0,0)),
	frequency_(1),
	amplitude_(0),
	phaseShift_(0)
{}


directionVaryingVelocityFvPatchVectorField::directionVaryingVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    meanValue_(dict.lookup("meanValue")),
	direction_(dict.lookup("direction")),
	frequency_(readScalar(dict.lookup("frequency"))),
	amplitude_(readScalar(dict.lookup("amplitude"))),
	phaseShift_(readScalar(dict.lookup("phaseShift")))
{
	direction_ /= mag(direction_);
	amplitude_ *= mathematicalConstant::pi/180.0;
	phaseShift_ *= mathematicalConstant::pi/180.0;

    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==(vectorField("value", dict, p.size()));
    }
    else
    {
        fixedValueFvPatchVectorField::operator==(vectorField(p.size(),currentValue()));
    }
}

directionVaryingVelocityFvPatchVectorField::directionVaryingVelocityFvPatchVectorField
(
    const directionVaryingVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    meanValue_(ptf.meanValue_),
	direction_(ptf.direction_),
	frequency_(ptf.frequency_),
	amplitude_(ptf.amplitude_),
	phaseShift_(ptf.phaseShift_)
{}

directionVaryingVelocityFvPatchVectorField::directionVaryingVelocityFvPatchVectorField
(
    const directionVaryingVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    meanValue_(ptf.meanValue_),
	direction_(ptf.direction_),
	frequency_(ptf.frequency_),
	amplitude_(ptf.amplitude_),
	phaseShift_(ptf.phaseShift_)
{}

directionVaryingVelocityFvPatchVectorField::directionVaryingVelocityFvPatchVectorField
(
    const directionVaryingVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    meanValue_(ptf.meanValue_),
	direction_(ptf.direction_),
	frequency_(ptf.frequency_),
	amplitude_(ptf.amplitude_),
	phaseShift_(ptf.phaseShift_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void directionVaryingVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}

void directionVaryingVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}

void directionVaryingVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

	vectorField& patchField = *this;
    patchField = vectorField(patchField.size(),currentValue());

    fixedValueFvPatchVectorField::updateCoeffs();
}

void directionVaryingVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("meanValue") << meanValue_ << token::END_STATEMENT << nl;
	os.writeKeyword("direction") << direction_ << token::END_STATEMENT << nl;
	os.writeKeyword("frequency") << frequency_ << token::END_STATEMENT << nl;
	os.writeKeyword("amplitude") << amplitude_*180.0/mathematicalConstant::pi << token::END_STATEMENT << nl;
	os.writeKeyword("phaseShift") << phaseShift_*180.0/mathematicalConstant::pi << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    directionVaryingVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
