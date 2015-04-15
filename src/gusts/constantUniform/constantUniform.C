
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
    constantUniform

\*---------------------------------------------------------------------------*/
#include "constantUniform.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(constantUniform, 0);
addToRunTimeSelectionTable(gust, constantUniform, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constantUniform::constantUniform(const dictionary& dict,const fvMesh& mesh):
gust(mesh),
dict_(dict),
mesh_(mesh),
startTime_(readScalar(dict.lookup("startTime"))),
endTime_(readScalar(dict.lookup("endTime"))),
gustAmplitude_(dict.lookup("amplitude")),
gustActive_(false),
smoothStart_(readBool(dict.lookup("smoothStart"))),
dT_(0)
{
    if(smoothStart_==true){
        Info << "Using smoothStart" << endl;
        dT_ = readScalar(dict.lookup("dT"));
    }
    Info << "Gust properties: T0 = " << startTime_ << ", T1 = " << endTime_ << ", dT = " << dT_ <<", A = " << gustAmplitude_.value() << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constantUniform::~constantUniform()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceVectorField& constantUniform::getGustVelocity(const scalar t,const volVectorField& U)
{
    dimensionedVector currentGustAmplitude = gustAmplitude_;
    if(t >= startTime_ && t <= endTime_)
    {
        gustActive_ = true;

        if(smoothStart_ && t<=(dT_+startTime_)){
            currentGustAmplitude.value() = gustAmplitude_.value()*0.5*(1.0 - Foam::cos(mathematicalConstant::pi/dT_*(t-startTime_)));
        }else{
            currentGustAmplitude.value() = gustAmplitude_.value();   
        }
    }else{
        gustActive_ = false;
        currentGustAmplitude.value() = vector::zero;
    }

    Ugust_ == currentGustAmplitude;

    Info << "Gust amplitude = " << currentGustAmplitude.value() << endl;

    return Ugust_;
}

const volVectorField& constantUniform::getGustVelocityCells(const scalar t,const volVectorField& U)
{
    dimensionedVector currentGustAmplitude = gustAmplitude_;
    if(t >= startTime_ && t <= endTime_)
    {
        gustActive_ = true;

        if(smoothStart_ && t<=(dT_+startTime_)){
            currentGustAmplitude.value() = gustAmplitude_.value()*0.5*(1.0 - Foam::cos(mathematicalConstant::pi/dT_*(t-startTime_)));
        }else{
            currentGustAmplitude.value() = gustAmplitude_.value();   
        }
    }else{
        gustActive_ = false;
        currentGustAmplitude.value() = vector::zero;
    }

    UgustCell_ == currentGustAmplitude;

    return UgustCell_;
}

const bool constantUniform::isActive() const
{
    return gustActive_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
