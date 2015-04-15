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
    twoDoFRigidBody

\*---------------------------------------------------------------------------*/
#include "twoDoFRigidBody.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "simpleMatrix.H"
#include "tensor.H"

using namespace Foam;

defineTypeNameAndDebug(twoDoFRigidBody, 0);
addToRunTimeSelectionTable(rigidBody, twoDoFRigidBody, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //
void twoDoFRigidBody::writeState(const scalar& time)
{
    //Body
    scalar bodyForce = bodyForces()[1];
    scalar bodyMoment = bodyMoments()[2];
    outputFile()    << time << "\t" << bodyState()[0] << "\t" << bodyState()[1] << "\t"
                    << bodyState()[2] << "\t" << bodyState()[3] << "\t"
                    << bodyForce << "\t" << bodyMoment << endl;
}

void twoDoFRigidBody::writeHeader()
{
    if(writeToFile())
    {
        outputFile() << "time\ty\ttheta\tdydt\tdthetadt\tFy\tMtheta" << endl;
    }
}


scalar twoDoFRigidBody::getForce(const vector forces) const
{
    return forces.y();
}

scalar twoDoFRigidBody::getMoment(const vector moments) const
{
    return moments.z();
}


void twoDoFRigidBody::calculateMoments(const Field<vectorField>& patchForces){
    scalar ybody = bodyState()[0];
    vector CofR = rc_ + vector(0.0,ybody,0.0);

    vector totalMoments(vector::zero);
    forAll(getPatchIDs(),ipatch)
    {
        vectorField Md = mesh().boundaryMesh()[getPatchIDs()[ipatch]].faceCentres() - CofR;
        vectorField moments = (Md ^ patchForces[getPatchIDs()[ipatch]]);

        //Parallel implementation
        totalMoments += gSum(moments);
    }

    bodyMoments() = totalMoments;
}

scalar twoDoFRigidBody::calculateThetaDamping(const scalar& theta) const{

    scalar ctheta = ctheta_;
    //Some function to calculate ctheta
    //ctheta = ctheta_*(1 + theta + pow(theta,2));

    return ctheta;
}


scalar twoDoFRigidBody::calculatePosition()
{
    Info << "Solve Solid" << endl;
    const scalar& dt = mesh().time().deltaT().value();

    //Start by saving previous iteration
    bodyStatePrevIter() = bodyState();

    //Set values of x,y,theta,dxdt,dydt,dthetadt to old values initially (every subiteration of fsi also)
    scalar ybody = bodyStateOld()[0];
    scalar theta = bodyStateOld()[1];
    scalar dybody = bodyStateOld()[2];
    scalar dtheta = bodyStateOld()[3];

    if(integrationScheme()=="RK4")
    {
        SquareMatrix<scalar> rkMatrix = getRK4Coeffs();

        label nStages=rkMatrix.n()-1;
        Field<scalarField> bodyStateVelStages(rkMatrix.n()-1,scalarField(bodyState().size(),0.0));
        for(int i=0;i<nStages;i++)
        {
            //Determine current state based on bodyStateVel
            scalarField bodyStageState(bodyStateOld());
            for(int k=0;k<i;k++){
                bodyStageState += dt*rkMatrix[i][k+1]*bodyStateVelStages[k];
            }
            ybody=bodyStageState[0];dybody=bodyStageState[2];
            theta=bodyStageState[1];dtheta=bodyStageState[3];

            //Interpolate forces
            vector bodyForces = interpolateRK4forces(rkMatrix[i][0]);
            scalar bodyForce = getForce(bodyForces);
            vector bodyMoments = interpolateRK4moments(rkMatrix[i][0]);
            scalar bodyMoment = getMoment(bodyMoments);

            //Get damping coefficient
            scalar ctheta = calculateThetaDamping(theta);

            //Calculate the acceleration of state variables
            simpleMatrix <scalar> C(2,0,0);
            C[0][0] = m_;
            C[0][1] = m_*Lcg2rc_*cos(theta+thetaGeom_);
            C[1][0] = m_*Lcg2rc_*cos(theta+thetaGeom_);
            C[1][1] = Icg_ + m_*pow(Lcg2rc_,2);

            C.source()[0] = bodyForce - c_*dybody - k_*ybody + m_*Lcg2rc_*pow(dtheta,2)*sin(theta+thetaGeom_);
            C.source()[1] = bodyMoment - ctheta*dtheta - ktheta_*theta;
            List<scalar> bodyStateAcc = C.LUsolve();

            //Calculate bodyStateVel for this stage.
            bodyStateVelStages[i][0]=dybody;
            bodyStateVelStages[i][1]=dtheta;
            bodyStateVelStages[i][2]=bodyStateAcc[0];
            bodyStateVelStages[i][3]=bodyStateAcc[1];
        }
        //calculate final state at t(n+1)
        bodyState() = bodyStateOld();
        for(int i=0;i<nStages;i++){
            bodyState() += dt*bodyStateVelStages[i]*rkMatrix[4][i+1];
        }
    }
    /*else if(integrationScheme()=="BDF1")
    {
        //Set old body state variables
        scalar ybody0 = bodyStateOld()[0];
        scalar dybody0 = bodyStateOld()[1];

        //Matrix with equations
        simpleMatrix<scalar> A(2,0,0);
        A[0][0] = c_/m_;
        A[0][1] = k_/m_;
        A[1][0] = -1;
        A[1][1] = 0;

        //Matrix with diagonal BDF1 terms
        A[0][0] += 1/dt;
        A[1][1] += 1/dt;

        scalar bodyForce = getForce(bodyForces());

        //rhs
        A.source()[0] = bodyForce/m_ - (-1/dt*dybody0);
        A.source()[1] = - (-1/dt*ybody0);

        //solve, where solution is dybody,ybody
        List<scalar> bodyState = A.LUsolve();
        this->bodyState()[1] = bodyState[0];
        this->bodyState()[0] = bodyState[1];
    }
    else if(integrationScheme()=="BDF2")
    {
        scalar deltaT = deltaT_();
        scalar deltaT0 = deltaT0_();
        scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
        scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
        scalar coefft0  = -1*(coefft + coefft00);

        //Set old body state variables
        scalar ybody0 = bodyStateOld()[0];
        scalar dybody0 = bodyStateOld()[1];
        scalar ybody00 = bodyStateOldOld()[0];
        scalar dybody00 = bodyStateOldOld()[1];

        //Matrix with equations
        simpleMatrix<scalar> A(2,0,0);
        A[0][0] = c_/m_;
        A[0][1] = k_/m_;
        A[1][0] = -1;
        A[1][1] = 0;

        //Matrix with diagonal BDF2 terms
        A[0][0] += coefft/dt;
        A[1][1] += coefft/dt;

        scalar bodyForce = getForce(bodyForces());

        //rhs
        A.source()[0] = bodyForce/m_ - (coefft0/dt*dybody0 + coefft00*dybody00/dt);
        A.source()[1] = - (coefft0/dt*ybody0 + coefft00*ybody00/dt);

        //solve, where solution is dybody,ybody
        List<scalar> bodyState = A.LUsolve();
        this->bodyState()[1] = bodyState[0];
        this->bodyState()[0] = bodyState[1];
    }
    else if(integrationScheme()=="BDF3")
    {
        scalarField coeffs = bdf3Coeffs();

        //Set old body state variables
        scalar ybody0 = bodyStateOld()[0];
        scalar dybody0 = bodyStateOld()[1];
        scalar ybody00 = bodyStateOldOld()[0];
        scalar dybody00 = bodyStateOldOld()[1];
        scalar ybody000 = bodyStateOldOldOld()[0];
        scalar dybody000 = bodyStateOldOldOld()[1];

        //Matrix with equations
        simpleMatrix<scalar> A(2,0,0);
        A[0][0] = c_/m_;
        A[0][1] = k_/m_;
        A[1][0] = -1;
        A[1][1] = 0;

        //Matrix with diagonal BDF2 terms
        A[0][0] += coeffs[0]/dt;
        A[1][1] += coeffs[0]/dt;

        scalar bodyForce = getForce(bodyForces());

        //rhs
        A.source()[0] = bodyForce/m_ - (coeffs[1]/dt*dybody0 + coeffs[2]*dybody00/dt + coeffs[3]*dybody000/dt);
        A.source()[1] = - (coeffs[1]/dt*ybody0 + coeffs[2]*ybody00/dt + coeffs[3]*ybody000/dt);

        //solve, where solution is dybody,ybody
        List<scalar> bodyState = A.LUsolve();
        this->bodyState()[1] = bodyState[0];
        this->bodyState()[0] = bodyState[1];
    }*/
    else
    {
        FatalErrorIn("void twoDoFRigidBody::calculatePosition()")
            << "integrationScheme " << integrationScheme() << " not available."
            << "Available schemes are: (RK4,BDF1, BDF2, BDF3)"
            << abort(FatalError);
    }

    return bodyState()[0];
}

void twoDoFRigidBody::readInitialStateIfPresent(const dictionary& dict)
{
    bodyState() = scalarField(4,0.0);

    scalar stateY = dict.lookupOrDefault("initialY",0.0);
    scalar stateTheta = dict.lookupOrDefault("initialTheta",0.0);
    scalar stateVelY = dict.lookupOrDefault("initialYVelocity",0.0);
    scalar stateVelTheta = dict.lookupOrDefault("initialThetaVelocity",0.0);

    bodyState()[0] = stateY;
    bodyState()[1] = stateTheta;
    bodyState()[2] = stateVelY;
    bodyState()[3] = stateVelTheta;

    bodyStateOld() = bodyState();
    bodyStateOldOld() = bodyState();
    bodyStateOldOldOld() = bodyState();
    bodyStatePrevIter() = bodyState();
    bodyState0() = bodyState();

    Info << "Initial state is set to: " << bodyState() << endl;
}


// * * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoDoFRigidBody::twoDoFRigidBody(const dictionary& dict,const fvMesh& mesh,const word name):
rigidBody(dict.parent(),mesh,name),
k_(readScalar(dict.lookup("k"))),
c_(readScalar(dict.lookup("c"))),
ktheta_(readScalar(dict.lookup("ktheta"))),
ctheta_(readScalar(dict.lookup("ctheta"))),
m_(0.0),
Icg_(readScalar(dict.lookup("Icg"))),
Lcg2rc_(readScalar(dict.lookup("Lcg2rc"))),
rc_(dict.lookup("rc")),
thetaGeom_(dict.lookupOrDefault("thetaGeom",0.0))
{
    //Determine thickness of mesh and set m accordingly
    Vector<label> e = (-mesh.geometricD()+Vector<label>(1,1,1))/2;
    scalar meshDepth(vector(e.x(),e.y(),e.z()) & vector(mesh.bounds().span()));
    scalar m(readScalar(dict.lookup("m")));
    m_ = m*meshDepth;

    Info << "BodyProperties: m = " << m_ << endl;

    //Read initial state
    readInitialStateIfPresent(dict);

    //Create header of output file in case writing is requested
    writeHeader();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

twoDoFRigidBody::~twoDoFRigidBody()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::septernion twoDoFRigidBody::calculateTransformation(const Field<vectorField>& patchForces)
{
    //- Update forces
    calculateForces(patchForces);
    calculateMoments(patchForces);

    //- Calculate new position based on new forces
    calculatePosition();

    vector translationVector(0.0,bodyState()[0]-bodyState0()[0],0.0);
    quaternion R(0,0,bodyState()[1]-bodyState0()[1]);
    septernion TR(septernion(rc_ + translationVector)*R*septernion(-rc_));

    return TR;
}

void twoDoFRigidBody::write()
{
    if(writeToFile())
    {
        writeState(mesh().time().value());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
