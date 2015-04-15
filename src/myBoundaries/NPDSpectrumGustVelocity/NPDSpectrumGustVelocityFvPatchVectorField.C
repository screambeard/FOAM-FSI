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

#include "NPDSpectrumGustVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mathematicalConstants.H"
#include "Random.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar NPDSpectrumGustVelocityFvPatchVectorField::currentAmplitude() const
{
	const scalar& pi = mathematicalConstant::pi;
    const scalar& currentTime = this->db().time().value();
    scalar currentAmplitude = 0.0;
    forAll(frequencies_,i){
        currentAmplitude += amplitudes_[i]*sin(2*pi*frequencies_[i]*currentTime + phaseShifts_[i]);
    }

    Info << "currentAmplitude = " << currentAmplitude << endl;
    return currentAmplitude;
}

void NPDSpectrumGustVelocityFvPatchVectorField::setSpectrum()
{
    if(Nsegments_>1){
        
        //Set frequencies with uniform distribution between max and min
        if(maxFrequency_<=minFrequency_){
            FatalErrorIn
            (
                "NPDSpectrumGustVelocityFvPatchVectorField::setSpectrum()"
            )   << "maxFrequency > minFrequency"
                << abort(FatalError);
        }
        
        scalar dF = (maxFrequency_ - minFrequency_)/(Nsegments_-1);
        forAll(frequencies_,i){
            frequencies_[i] = minFrequency_ + dF*i;
        }

        //Amplitudes (based on Uref = 10 m/s)
        forAll(amplitudes_,i){
            amplitudes_[i] = getAmplitude(frequencies_[i]);
        }

        //PhaseShifts (random generated between 0 and 2pi)
        forAll(phaseShifts_,i){
            phaseShifts_[i] = rndGen_.scalar01()*2*mathematicalConstant::pi;
            Info << "phaseShift segment" << i << " = " << phaseShifts_[i] << "(" << phaseShifts_[i]*180/mathematicalConstant::pi << ")" <<  endl;
        }

    }else if(Nsegments_ ==1){

    }else{
        FatalErrorIn
        (
            "NPDSpectrumGustVelocityFvPatchVectorField::setSpectrum()"
        )   << "Nsegments must be bigger than 0"
            << abort(FatalError);
    }
}

scalar NPDSpectrumGustVelocityFvPatchVectorField::getAmplitude(scalar frequency) const
{
    scalar n=0.468;
    scalar fbar=172*frequency*pow(Uref_/10,-3.0/4.0);
    scalar amplitude = (320*pow(Uref_/10,2))/(pow(1+pow(fbar,n),5/3/n));
    return amplitude;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NPDSpectrumGustVelocityFvPatchVectorField::NPDSpectrumGustVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_(p.size()),
    amplitude_(p.size()),
    direction_(vector::zero),
    minFrequency_(0.0),
    maxFrequency_(0.0),
    Nsegments_(0),
    Uref_(10),
    phaseShifts_(Nsegments_,0.0),
    frequencies_(Nsegments_,0.0),
    amplitudes_(Nsegments_,0.0),
    rndGen_(1),
    curTimeIndex_(-1)
{}


NPDSpectrumGustVelocityFvPatchVectorField::NPDSpectrumGustVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    refValue_("refValue", dict, p.size()),
    amplitude_("amplitude",dict, p.size()),
    direction_(dict.lookup("direction")),
    minFrequency_(readScalar(dict.lookup("minFrequency"))),
    maxFrequency_(readScalar(dict.lookup("maxFrequency"))),
    Nsegments_(readLabel(dict.lookup("Nsegments"))),
    Uref_(readScalar(dict.lookup("Uref"))),
    phaseShifts_(Nsegments_,0.0),
    frequencies_(Nsegments_,0.0),
    amplitudes_(Nsegments_,0.0),
    rndGen_(1),
    curTimeIndex_(-1)
{
    setSpectrum();

    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==(vectorField("value", dict, p.size()));
    }
    else
    {
        fixedValueFvPatchVectorField::operator==(refValue_+direction_*currentAmplitude());
    }
}

NPDSpectrumGustVelocityFvPatchVectorField::NPDSpectrumGustVelocityFvPatchVectorField
(
    const NPDSpectrumGustVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    refValue_(ptf.refValue_, mapper),
    amplitude_(ptf.amplitude_, mapper),
    direction_(ptf.direction_),
    minFrequency_(ptf.minFrequency_),
    maxFrequency_(ptf.maxFrequency_),
    Nsegments_(ptf.Nsegments_),
    Uref_(ptf.Uref_),
    phaseShifts_(ptf.phaseShifts_),
    frequencies_(ptf.frequencies_),
    amplitudes_(ptf.amplitudes_),
    rndGen_(ptf.rndGen_),
    curTimeIndex_(-1)
{}

NPDSpectrumGustVelocityFvPatchVectorField::NPDSpectrumGustVelocityFvPatchVectorField
(
    const NPDSpectrumGustVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    refValue_(ptf.refValue_),
    amplitude_(ptf.amplitude_),
    direction_(ptf.direction_),
    minFrequency_(ptf.minFrequency_),
    maxFrequency_(ptf.maxFrequency_),
    Nsegments_(ptf.Nsegments_),
    Uref_(ptf.Uref_),
    phaseShifts_(ptf.phaseShifts_),
    frequencies_(ptf.frequencies_),
    amplitudes_(ptf.amplitudes_),
    rndGen_(ptf.rndGen_),
    curTimeIndex_(-1)
{}

NPDSpectrumGustVelocityFvPatchVectorField::NPDSpectrumGustVelocityFvPatchVectorField
(
    const NPDSpectrumGustVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    refValue_(ptf.refValue_),
    amplitude_(ptf.amplitude_),
    direction_(ptf.direction_),
    minFrequency_(ptf.minFrequency_),
    maxFrequency_(ptf.maxFrequency_),
    Nsegments_(ptf.Nsegments_),
    Uref_(ptf.Uref_),
    phaseShifts_(ptf.phaseShifts_),
    frequencies_(ptf.frequencies_),
    amplitudes_(ptf.amplitudes_),
    rndGen_(ptf.rndGen_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void NPDSpectrumGustVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    refValue_.autoMap(m);
    amplitude_.autoMap(m);
}

void NPDSpectrumGustVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const NPDSpectrumGustVelocityFvPatchVectorField& tiptf =
        refCast<const NPDSpectrumGustVelocityFvPatchVectorField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
    amplitude_.rmap(tiptf.amplitude_, addr);
}

void NPDSpectrumGustVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        vectorField& patchField = *this;

        patchField = refValue_+direction_*currentAmplitude();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}

void NPDSpectrumGustVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    refValue_.writeEntry("refValue", os);
    os.writeKeyword("maxFrequency")
        << maxFrequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("minFrequency")
        << minFrequency_ << token::END_STATEMENT << nl;
    os.writeKeyword("Nsegments")
        << Nsegments_ << token::END_STATEMENT << nl;
    os.writeKeyword("Uref")
        << Uref_ << token::END_STATEMENT << nl;
    os.writeKeyword("direction")
        << direction_ << token::END_STATEMENT << nl;        
    os.writeKeyword("phaseShifts")
        << phaseShifts_ << token::END_STATEMENT << nl;
    os.writeKeyword("amplitudes")
        << amplitudes_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    NPDSpectrumGustVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
