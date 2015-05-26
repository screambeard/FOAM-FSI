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

#include "constantAcceleration.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(constantAcceleration, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        constantAcceleration,
        dictionary
    );
};
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector
Foam::solidBodyMotionFunctions::constantAcceleration::calcPosition
(
    const scalar t
) const
{
    vector position = 0.5*acceleration_*pow(t,2);
    Info << "positions = " << position << endl;

    return position;
}

Foam::vector
Foam::solidBodyMotionFunctions::constantAcceleration::calcVelocity
(
    const scalar t
) const
{
    vector velocity = acceleration_*t;

    return velocity;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::constantAcceleration::constantAcceleration
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    acceleration_(SBMFCoeffs_.lookup("acceleration"))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::constantAcceleration::~constantAcceleration()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::constantAcceleration::transformation() const
{
    scalar t = time_.value();

    septernion TR(calcPosition(t), quaternion::I);

    Info<< "solidBodyMotionFunctions::constantAcceleration::transformation(): "
        << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


Foam::septernion
Foam::solidBodyMotionFunctions::constantAcceleration::velocity() const
{
    scalar t = time_.value();
    scalar dt = time_.deltaT().value();

    septernion TV
    (
        (calcPosition(t + dt) - calcPosition(t))/dt,
        //calcVelocity(t),
        quaternion::zero
    );

    return TV;
}


bool Foam::solidBodyMotionFunctions::constantAcceleration::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    SBMFCoeffs_.lookup("acceleration") >> acceleration_;

    return true;
}


// ************************************************************************* //
