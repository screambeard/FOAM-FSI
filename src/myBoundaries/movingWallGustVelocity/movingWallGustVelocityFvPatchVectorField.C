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

#include "movingWallGustVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

movingWallGustVelocityFvPatchVectorField::movingWallGustVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero)
{}


movingWallGustVelocityFvPatchVectorField::movingWallGustVelocityFvPatchVectorField
(
    const movingWallGustVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    myTimeIndex_(ptf.myTimeIndex_),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero)
{}


movingWallGustVelocityFvPatchVectorField::movingWallGustVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    myTimeIndex_(dimensionedInternalField().mesh().time().timeIndex()),
    Fc_(p.patch().size(),vector::zero),
    oldFc_(p.patch().size(),vector::zero),
    oldoldFc_(p.patch().size(),vector::zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    Fc_ = p.patch().faceCentres();
    oldFc_ = p.patch().faceCentres();
    oldoldFc_ = p.patch().faceCentres();
}


movingWallGustVelocityFvPatchVectorField::movingWallGustVelocityFvPatchVectorField
(
    const movingWallGustVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.Fc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_)
{}


movingWallGustVelocityFvPatchVectorField::movingWallGustVelocityFvPatchVectorField
(
    const movingWallGustVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    myTimeIndex_(pivpvf.myTimeIndex_),
    Fc_(pivpvf.oldFc_),
    oldFc_(pivpvf.oldFc_),
    oldoldFc_(pivpvf.oldoldFc_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void movingWallGustVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //Get current mesh data
    const fvMesh& mesh = dimensionedInternalField().mesh();
    const fvPatch& p = patch();
    const polyPatch& pp = p.patch();
    
    //Update the face center locations in time
    if(myTimeIndex_< mesh.time().timeIndex())
    { 
        oldoldFc_ = oldFc_;
        oldFc_ = Fc_;
        Fc_ = pp.faceCentres();
        myTimeIndex_ = mesh.time().timeIndex();
    }

    const volVectorField& U = mesh.lookupObject<volVectorField>(dimensionedInternalField().name());
    const surfaceVectorField& Ugust = mesh.lookupObject<surfaceVectorField>(dimensionedInternalField().name()+"gust");
    const vectorField UgustPatch = p.patchField<surfaceVectorField, vector>(Ugust);

    scalar deltaT = mesh.time().deltaT().value();
    scalar deltaT0 = mesh.time().deltaT0().value();
    if (U.oldTime().timeIndex() == U.oldTime().oldTime().timeIndex() || U.oldTime().oldTime().timeIndex() < 0)
    {
        deltaT0 = GREAT;
    }

    //Set coefficients based on deltaT and deltaT0
    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    vectorField Up = (coefft*Fc_ - coefft0*oldFc_ + coefft00*oldoldFc_)/mesh.time().deltaT().value();
    
    scalarField phip = p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U));

    vectorField n = p.nf();
    const scalarField& magSf = p.magSf();
    scalarField Un = phip/(magSf + VSMALL);

    //Thijs Gillebaart 21/2/2014: not sure whether UgustPatch must be on both parallel and normal velocity or only in Un
    //vectorField::operator=(Up + n*(Un - (n & Up)));
    vectorField::operator=(-UgustPatch + Up + n*(Un - (n & Up)));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void movingWallGustVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    movingWallGustVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
