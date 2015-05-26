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

Class
    Foam::solidBodyMotionFunctions::constantAcceleration

Description
    Linear oscillation motion function

SourceFiles
    constantAcceleration.C

\*---------------------------------------------------------------------------*/

#ifndef constantAcceleration_H
#define constantAcceleration_H

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                      Class constantAcceleration Declaration
\*---------------------------------------------------------------------------*/

class constantAcceleration
:
    public solidBodyMotionFunction
{
    // Private data

        //- Motion amplitude
        vector acceleration_;

    // Private Member Functions

        //- Disallow copy construct
        constantAcceleration(const constantAcceleration&);

        //- Disallow default bitwise assignment
        void operator=(const constantAcceleration&);


        //- Calculate position
        vector calcPosition(const scalar t) const;

        //- Calculate position
        vector calcVelocity(const scalar t) const;


public:

    //- Runtime type information
    TypeName("constantAcceleration");


    // Constructors

        //- Construct from components
        constantAcceleration
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Destructor

        virtual ~constantAcceleration();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Return the solid-body motion velocity
        virtual septernion velocity() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
