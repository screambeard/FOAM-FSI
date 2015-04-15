
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    cosineUniform

\*---------------------------------------------------------------------------*/
#include "cosineUniform.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(cosineUniform, 0);
addToRunTimeSelectionTable(gust, cosineUniform, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cosineUniform::cosineUniform(const dictionary& dict,const fvMesh& mesh):
gust(mesh),
dict_(dict),
mesh_(mesh),
startTime_(readScalar(dict.lookup("startTime"))),
frequency_(readScalar(dict.lookup("frequency"))),
gustAmplitude_(dict.lookup("amplitude")),
gustActive_(false)
{
    if(frequency_<=0.0){
        FatalErrorIn("cosineUniform::cosineUniform(const dictionary& dict,const fvMesh& mesh)")
            << "frequency must be bigger than zero, instead of " << frequency_
            << abort(FatalError);
    }
    Info << "Gust properties: T0 = " << startTime_ << ", f = " << frequency_ <<", A = " << gustAmplitude_.value() << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cosineUniform::~cosineUniform()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceVectorField& cosineUniform::getGustVelocity(const scalar t,const volVectorField& U)
{
    dimensionedVector currentGustAmplitude = gustAmplitude_;
    if(t >= startTime_ && t <= startTime_+1.0/frequency_)
    {
        gustActive_ = true;
        currentGustAmplitude.value() = gustAmplitude_.value()*0.5*(1.0 - Foam::cos(2*mathematicalConstant::pi * frequency_ * (t-startTime_)));
    }else{
        gustActive_ = false;
        currentGustAmplitude.value() = vector::zero;
    }

    Ugust_ == currentGustAmplitude;

    Info << "Gust amplitude = " << currentGustAmplitude.value() << endl;

    return Ugust_;
}

const volVectorField& cosineUniform::getGustVelocityCells(const scalar t,const volVectorField& U)
{
    dimensionedVector currentGustAmplitude = gustAmplitude_;
    if(t >= startTime_ && t <= startTime_+1.0/frequency_)
    {
        gustActive_ = true;
        currentGustAmplitude.value() = gustAmplitude_.value()*0.5*(1.0 - Foam::cos(2*mathematicalConstant::pi * frequency_ * (t-startTime_)));
    }else{
        gustActive_ = false;
        currentGustAmplitude.value() = vector::zero;
    }

    UgustCell_ == currentGustAmplitude;

    return UgustCell_;
}

const bool cosineUniform::isActive() const
{
    return gustActive_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
