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
    gust

Description
    Body motion virtual base class

Author
    Thijs Gillebaart, TU Delft.  All rights reserved.

SourceFiles
    gust.C
    newgust.C

\*---------------------------------------------------------------------------*/

#ifndef gust_H
#define gust_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "primitiveFields.H"
#include "fvCFD.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gust Declaration
\*---------------------------------------------------------------------------*/

class gust
{
    // Private Member Functions

        //- Disallow copy construct
        gust(const gust&);

        //- Disallow default bitwise assignment
        void operator=(const gust&);

protected:

        //Velocity field gust
        surfaceVectorField Ugust_;

        //Velocity field gust at cell centres
        volVectorField UgustCell_;

public:

    //- Runtime type information
    TypeName("gust");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            gust,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Selectors

        //- Return a pointer to the selected RigidMotion function
        static autoPtr<gust> New
        (
            const word& type,
            const dictionary& dict,
			const fvMesh& mesh
        );


    // Constructors

        gust(const fvMesh& mesh);

        //- Create and return a clone
        virtual autoPtr<gust> clone() const = 0;


    // Destructor

        virtual ~gust()
        {}


    // Member Functions

        //- Returns gust velocity field
        virtual const surfaceVectorField& getGustVelocity(const scalar t,const volVectorField& U) = 0;

        //- Returns gust velocity field at cell centres
        virtual const volVectorField& getGustVelocityCells(const scalar t,const volVectorField& U) = 0;

        //- Return if gust is active
        virtual const bool isActive() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
