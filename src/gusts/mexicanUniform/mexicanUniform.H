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
    mexicanUniform

Description
	Rigid body motion based sinus functions

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    mexicanUniform.C

\*---------------------------------------------------------------------------*/

#ifndef mexicanUniform_H
#define mexicanUniform_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "gust.H"
#include "dimensionedVector.H"
#include "dimensionedScalar.H"
#include "Time.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mexicanUniform Declaration
\*---------------------------------------------------------------------------*/

class mexicanUniform
:
	public gust
{
		// Private data
		const dictionary& dict_;
		const fvMesh& mesh_;

        scalar startTime_;
        scalar frequency_;
        dimensionedVector gustAmplitude_;
        bool gustActive_;

        //- Writing function
        void writePosition(const scalar& time);

    	//- Disallow default bitwise copy construct
        mexicanUniform(const mexicanUniform&);
	    
	public:

    //- Runtime type information
    TypeName("mexicanUniform");

	//Public data

    // Constructors

        //- Construct from dict + mesh
		mexicanUniform(const dictionary& dict,const fvMesh& mesh);

		virtual autoPtr<gust> clone() const
        {
            return autoPtr<gust>(new mexicanUniform(this->dict_,this->mesh_));
        };

    // Destructor

        virtual ~mexicanUniform();

    // Member Functions

        //- Returns gust velocity field
        virtual const surfaceVectorField& getGustVelocity(const scalar t,const volVectorField& U);

        //- Returns gust velocity field at cell centres
        virtual const volVectorField& getGustVelocityCells(const scalar t,const volVectorField& U);

        //- Return if gust is active
        virtual const bool isActive() const;
};
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
