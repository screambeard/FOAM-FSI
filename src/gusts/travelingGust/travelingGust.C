
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
    travelingGust

\*---------------------------------------------------------------------------*/
#include "travelingGust.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(travelingGust, 0);
addToRunTimeSelectionTable(gust, travelingGust, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
const dimensionedVector travelingGust::getLocalAmplitude(const scalar relX) const{
    dimensionedVector currentGustAmplitude = gustAmplitude_;
    if(relX >= 0 && relX <= gustLength_){
    	if(gusttype_=="cosine"){
        	currentGustAmplitude.value() = gustAmplitude_.value()*0.5*(1-Foam::cos(2*mathematicalConstant::pi/gustLength_*relX));
    	}else if(gusttype_=="mexican"){
        	currentGustAmplitude.value() = -gustAmplitude_.value()*0.5*(1-Foam::cos(2*mathematicalConstant::pi/gustLength_*relX))*Foam::sin(3*mathematicalConstant::pi/gustLength_*relX);
        }else{
        	FatalErrorIn("dimensionedVector travelingGust::getLocalAmplitude(const scalar relX)")
		    << "Gusttype " << gusttype_ << " is not an option. Options are (cosine, mexican)"
		    << abort(FatalError);
        }
    }else{
        currentGustAmplitude.value() = gustAmplitude_.value()*0.0;
    }
    return currentGustAmplitude;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

travelingGust::travelingGust(const dictionary& dict,const fvMesh& mesh):
gust(mesh),
dict_(dict),
mesh_(mesh),
frequency_(readScalar(dict.lookup("frequency"))),
convectionSpeed_(dict.lookup("convectionSpeed")),
gustAmplitude_(dict.lookup("amplitude")),
origin_(dict.lookup("origin")),
gusttype_(dict.lookup("type")),
gustLength_(mag(convectionSpeed_)/frequency_),
gustActive_(true)
{
    Info << "Gust properties of travelingGust: f = " << frequency_ << ", Uconv = " << convectionSpeed_  <<", origin = " << origin_ <<", A = " << gustAmplitude_.value() << endl;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

travelingGust::~travelingGust()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceVectorField& travelingGust::getGustVelocity(const scalar t,const volVectorField& U)
{
    vector gustFront = convectionSpeed_*t+origin_;
    Info << "Gust front at " << gustFront << endl;

    //Do all internal faces
    forAll(Ugust_,iFace){
        const vector& position(mesh_.faceCentres()[iFace]);
        vector currentRelPos = position - gustFront + gustLength_*convectionSpeed_/mag(convectionSpeed_);
        //scalar currentRelX = position.x() - xGustFront + gustLength_;
        scalar currentRelX = (currentRelPos & (convectionSpeed_/(mag(convectionSpeed_))));

        Ugust_[iFace] = getLocalAmplitude(currentRelX).value();
    }

    //Do all boundary faces
    forAll(Ugust_.boundaryField(),iBoundary){
        const scalar& startInd = mesh_.boundaryMesh()[iBoundary].start();
        forAll(Ugust_.boundaryField()[iBoundary],iFace){
            const vector& position(mesh_.faceCentres()[iFace+startInd]);
            vector currentRelPos = position - gustFront + gustLength_*convectionSpeed_/mag(convectionSpeed_);
            //scalar currentRelX = position.x() - xGustFront + gustLength_;
            scalar currentRelX = (currentRelPos & (convectionSpeed_/(mag(convectionSpeed_))));

            Ugust_.boundaryField()[iBoundary][iFace] = getLocalAmplitude(currentRelX).value();
        }
    }

    return Ugust_;
}

const volVectorField& travelingGust::getGustVelocityCells(const scalar t,const volVectorField& U)
{
    vector gustFront = convectionSpeed_*t+origin_;
    Info << "Gust front at " << gustFront << endl;

    //Do all cells
    forAll(UgustCell_,iCell){
        const vector& position(mesh_.cellCentres()[iCell]);
        vector currentRelPos = position - gustFront + gustLength_*convectionSpeed_/mag(convectionSpeed_);
        
        scalar currentRelX = (currentRelPos & (convectionSpeed_/(mag(convectionSpeed_))));

        UgustCell_[iCell] = getLocalAmplitude(currentRelX).value();
    }

    //Do all boundary faces
    forAll(UgustCell_.boundaryField(),iBoundary){
        const scalar& startInd = mesh_.boundaryMesh()[iBoundary].start();
        forAll(UgustCell_.boundaryField()[iBoundary],iFace){
            const vector& position(mesh_.faceCentres()[iFace+startInd]);
            vector currentRelPos = position - gustFront + gustLength_*convectionSpeed_/mag(convectionSpeed_);

            scalar currentRelX = (currentRelPos & (convectionSpeed_/(mag(convectionSpeed_))));

            UgustCell_.boundaryField()[iBoundary][iFace] = getLocalAmplitude(currentRelX).value();
        }
    }

    return UgustCell_;
}

const bool travelingGust::isActive() const
{
    return gustActive_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
