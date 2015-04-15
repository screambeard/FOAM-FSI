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

\*---------------------------------------------------------------------------*/

#include "gust.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(gust, 0);
defineRunTimeSelectionTable(gust, dictionary);

gust::gust(const fvMesh& mesh):
Ugust_(IOobject("Ugust",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::AUTO_WRITE),mesh,dimensionSet(0,1,-1,0,0,0,0),fixedValueFvPatchVectorField::typeName),
UgustCell_(IOobject("UgustCell",mesh.time().timeName(),mesh,IOobject::NO_READ,IOobject::NO_WRITE),mesh,dimensionSet(0,1,-1,0,0,0,0),fixedValueFvPatchVectorField::typeName)
{

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
