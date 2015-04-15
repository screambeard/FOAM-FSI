
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
    travelingVortex

\*---------------------------------------------------------------------------*/
#include "travelingVortex.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(travelingVortex, 0);
addToRunTimeSelectionTable(gust, travelingVortex, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

travelingVortex::travelingVortex(const dictionary& dict,const fvMesh& mesh):
gust(mesh),
dict_(dict),
mesh_(mesh),
R_(readScalar(dict.lookup("radius"))),
strength_(readScalar(dict.lookup("strength"))),
convectionSpeed_(dict.lookup("convectionSpeed")),
origin_(dict.lookup("origin"))
{
    Info << "Gust properties of travelingVortex: R = " << R_ <<", Strength = " << strength_ <<", Uconv = " << convectionSpeed_  <<", origin = " << origin_ << ", Umax = " << Foam::exp(-0.5)*strength_/R_ <<  endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

travelingVortex::~travelingVortex()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceVectorField& travelingVortex::getGustVelocity(const scalar t,const volVectorField& U)
{
    vector gustOrigin = convectionSpeed_*t+origin_;
    Info << "Gust origin at " << gustOrigin << endl;

    //Calculate radius and all face values (internal & boundary)
    scalarField radius = mag(mesh_.faceCentres() - gustOrigin);
    vectorField UgustTmp(mesh_.nFaces(),vector::zero);
    UgustTmp.replace(0,-strength_*(mesh_.faceCentres().component(1) - gustOrigin.y())/sqr(R_) * Foam::exp(-0.5*sqr(radius/R_)));
    UgustTmp.replace(1,strength_*(mesh_.faceCentres().component(0) - gustOrigin.x())/sqr(R_) * Foam::exp(-0.5*sqr(radius/R_)));

    //Set internal field
    Ugust_.internalField() = SubList<vector>(UgustTmp,mesh_.nInternalFaces(),0);
    Info << "max(Ugust) = " << max(mag(Ugust_.internalField())) << endl;

    //Do all boundary faces
    forAll(Ugust_.boundaryField(),iBoundary){
        if(mesh_.boundaryMesh().types()[iBoundary] != "empty")
        {
            const scalar& startInd = mesh_.boundaryMesh()[iBoundary].start();
            const label& patchSize = mesh_.boundaryMesh()[iBoundary].size();
            Ugust_.boundaryField()[iBoundary] == vectorField(SubList<vector>(UgustTmp,patchSize,startInd));
            //Info << "Ugust at boundary " << iBoundary << ": " << Ugust_.boundaryField()[iBoundary] << endl;
        }
        
    }

    return Ugust_;
}

const volVectorField& travelingVortex::getGustVelocityCells(const scalar t,const volVectorField& U)
{
    return UgustCell_;
}

const bool travelingVortex::isActive() const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
