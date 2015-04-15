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
    rigidBody

\*---------------------------------------------------------------------------*/
#include "rigidBody.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(rigidBody, 0);
defineRunTimeSelectionTable(rigidBody, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
void rigidBody::checkInterpolationRK4force()
{
    if(integrationScheme_=="RK4"){
        switch(interpolationRK4force_)
        {
            case 1:
                Info << "rigidBody: Doing force interpolation at t(n+1)" << endl;
                break;
            case 2:
                Info << "rigidBody: Doing force interpolation at 0.5*[t(n+1)+t(n)]" << endl;
                break;
            case 4:
                Info << "rigidBody: Doing force interpolation at [k*t(n+1) + (1-k)*t(n)]" << endl;
                break;
            default:
                Info << "rigidBody: Unknown force interpolation for RK4. Setting to best." << endl;
                interpolationRK4force_ = 4;
                break;
        }
    }
}

void rigidBody::setPatchIDs(){
    forAll(patchNames_,ipatch){
        label patchID(mesh_.boundaryMesh().findPatchID(patchNames_[ipatch]));
        if(patchID < 0){
            FatalErrorIn("void rigidBody::setPatchIDs()")
            << "patchName " << patchNames_[ipatch] << " does not exist."
            << abort(FatalError);
        }else{
            patchIDs_[ipatch] = patchID;
        }
    }
}


vector rigidBody::interpolateRK4term(const scalar rkCoeff,const vector valueOld, const vector value) const
{
    vector termRK4 = vector::zero;
    switch(interpolationRK4force_)
    {
        case 1:
            termRK4 = value;
            break;
        case 2:
            termRK4 = 0.5*(value + valueOld);
            break;
        case 4:
            termRK4 = rkCoeff*value + (1-rkCoeff)*valueOld;
            break;
        default:
            termRK4 = rkCoeff*value + (1-rkCoeff)*valueOld;
            break;
        break;
    }

    return termRK4;
}

// * * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //

const vector& rigidBody::calculateForces(const Field<vectorField>& patchForces)
{
    vector forces = vector::zero;
    forAll(patchIDs_,ipatch)
    {
        forces += gSum(patchForces[patchIDs_[ipatch]]);
    }

    bodyForces_ = forces;

    return bodyForces_;
}

vector rigidBody::interpolateRK4forces(const scalar rkCoeff) const
{
    return interpolateRK4term(rkCoeff,bodyForcesOld_,bodyForces_);
}

vector rigidBody::interpolateRK4moments(const scalar rkCoeff) const
{
    return interpolateRK4term(rkCoeff,bodyMomentsOld_,bodyMoments_);
}

SquareMatrix<scalar> rigidBody::getRK4Coeffs() const
{
    // === matrix 5x5 === //
    // 0     0     0     0     0
    // 1/3   1/3   0     0     0
    // 2/3   -1/3  1     0     0
    // 1     1     -1    1     0
    // 0     1/8   3/8   3/8   1/8
    // === end === //
    SquareMatrix<scalar> rkMatrix(5,0.0);
    rkMatrix[1][0]=1.0/3.0;rkMatrix[2][0]=2.0/3.0;rkMatrix[3][0]=1.0;
    rkMatrix[4][1]=1.0/8.0;rkMatrix[4][2]=3.0/8.0;rkMatrix[4][3]=3.0/8.0;rkMatrix[4][4]=1.0/8.0;
    rkMatrix[1][1]=1.0/3.0;
    rkMatrix[2][1]=-1.0/3.0;rkMatrix[2][2]=1.0;
    rkMatrix[3][1]=1.0;rkMatrix[3][2]=-1.0;rkMatrix[3][3]=1.0;

    /*// === matrix 5x5 === //
    // 0     0     0     0     0
    // 1/2   1/2   0     0     0
    // 1/2   0     1/2   0     0
    // 1     0     0     1     0
    // 0     1/6   1/3   1/3   1/6
    // === end === //
    SquareMatrix<scalar> rkMatrix(5,0.0);
    rkMatrix[1][0]=1.0/2.0;rkMatrix[2][0]=1.0/2.0;rkMatrix[3][0]=1.0;
    rkMatrix[4][1]=1.0/6.0;rkMatrix[4][2]=1.0/3.0;rkMatrix[4][3]=1.0/3.0;rkMatrix[4][4]=1.0/6.0;
    rkMatrix[1][1]=1.0/2.0;
    rkMatrix[2][2]=1.0/2.0;
    rkMatrix[3][3]=1.0;*/

    return rkMatrix;
}

scalarField rigidBody::bdf3Coeffs() const
{
    scalarField coeffs(4,0.0);

    //if first time step
    if (mesh_.time().timeIndex() < 2)
    {
        coeffs[0]=1.0;
        coeffs[1]=-1.0;
    }
    else if(mesh_.time().timeIndex() < 3)
    {
        coeffs[0]=3.0/2.0;
        coeffs[1]=-2.0;
        coeffs[2]=1.0/2.0;
    }
    else
    {
        coeffs[0]=11.0/6.0;
        coeffs[1]=-3.0;
        coeffs[2]=3.0/2.0;
        coeffs[3]=-1.0/3.0;
    }

    return coeffs;
}

scalar rigidBody::deltaT_() const
{
    return mesh_.time().deltaT().value();
}


scalar rigidBody::deltaT0_() const
{
    if (mesh_.time().timeIndex() < 2)
    {
        Info << nl << "rigidBody: Doing Euler time integration in first time-step" << nl <<endl;
        return GREAT;
    }
    else
    {
        return mesh_.time().deltaT0().value();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rigidBody::rigidBody(const dictionary& dict,const fvMesh& mesh,const word name):
dict_(dict),
mesh_(mesh),
name_(name),
timeIndex_(-1),
write_(dict.lookupOrDefault("write",true)),
patchNames_(dict.lookup("patches")),
patchIDs_(patchNames_.size(),-1),
integrationScheme_(dict_.lookup("integrationScheme")),
interpolationRK4force_(dict_.lookupOrDefault("interpolation",4)),
bodyState_(0,0.0),
bodyStatePrevIter_(0,0.0),
bodyStateOld_(0,0.0),
bodyStateOldOld_(0,0.0),
bodyStateOldOldOld_(0,0.0),
bodyState0_(0,0.0),
bodyForces_(vector::zero),
bodyMoments_(vector::zero),
bodyForcesOld_(vector::zero),
bodyMomentsOld_(vector::zero),
ofBody_("rigidBody-"+name+"-state.dat")
{
    setPatchIDs();
    checkInterpolationRK4force();

    if(!writeToFile()){
        rm("body-"+name+"-state.dat");
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

rigidBody::~rigidBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& rigidBody::getPatchIDs() const{
	return patchIDs_;
}

void rigidBody::update()
{
    if(localTimeIndex() < mesh_.time().timeIndex()){
        bodyStateOldOldOld_ = bodyStateOldOld_;//set t(n-2) bodyState
        bodyStateOldOld_ = bodyStateOld_;//set t(n-1) bodyState
        bodyStateOld_ = bodyState_; //Set t(n) bodyState
        bodyForcesOld_ = bodyForces_;
        bodyMomentsOld_ = bodyMoments_;
    }
    timeIndex_ = mesh_.time().timeIndex();
}

bool rigidBody::writeToFile() const{
    return write_;
}


const label& rigidBody::localTimeIndex() const{
	return timeIndex_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
