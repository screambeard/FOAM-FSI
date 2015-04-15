
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
    vortexParticleGust

\*---------------------------------------------------------------------------*/
#include "vortexParticleGust.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

using namespace Foam;

defineTypeNameAndDebug(vortexParticleGust, 0);
addToRunTimeSelectionTable(gust, vortexParticleGust, dictionary);

// * * * * * * * * * * * * * * * * Private Member Functions* * * * * * * * * //

void vortexParticleGust::initializeParticles()
{
    Info << "Start initializeParticles() " << endl;

    const scalar& pi = mathematicalConstant::pi;
    
    //Distrubute particles according to specified radial and circular number of particles
    scalar deltaR = particlesDistributionRadius_/NvRadial_;
    scalar deltaC = 2*pi/NvCirc_;

    for(label ir = 0;ir<NvRadial_;ir++){
        scalar radius = (ir+1)*deltaR;
        for(label ic = 0;ic<NvCirc_;ic++){
            label ip = (ir*NvCirc_)+ic;
            scalar angle = ic*deltaC;
            particlePosition_[ip+1] = origin_ + radius*vector(cos(angle),sin(angle),0.0);
            scalar particleArea = deltaC*deltaR*radius;
            scalar particleVorticity = strength_/sqr(R_)*exp(-0.5*sqr(radius/R_))*(2-sqr(radius/R_));
            particleStrength_[ip+1] = particleVorticity * particleArea;
            particleSize_[ip+1] = 2*deltaR;
        }
    }

    //Do first central particle
    particlePosition_[0] = origin_;
    particleStrength_[0] = 2*strength_/sqr(R_) * (pi*sqr(0.5*deltaR));
    particleSize_[0] = 2*deltaR;
}

void vortexParticleGust::moveParticles(const volVectorField& U,const scalar dt)
{
    Info << "Start moveParticles() " << endl;
    const scalar& pi = mathematicalConstant::pi;

    //Apply RK4 time integration for all particles
    List<scalar> rkCoeffs(4,1.0);
    rkCoeffs[0]=1.0/4.0;rkCoeffs[1]=1.0/3.0;rkCoeffs[2]=1.0/2.0;

    const vectorField oldParticlePosition = particlePosition_;

    vectorField bla(Nv_,vector::zero);
    forAll(rkCoeffs,iStage)
    {
        //Interpolate "freestream" velocity from U @ particle positions
        vectorField particleVelocity(Nv_,vector(1.0,0.0,0.0));

        forAll(particlePosition_,iParticle)
        {
            const scalar& pstr = particleStrength_[iParticle];
            const vector& ppos = particlePosition_[iParticle];
            const scalar& psize = particleSize_[iParticle];
            const scalarField radius = mag(particlePosition_ - ppos)+SMALL;

            bla.replace(0,-(particlePosition_.component(1)-ppos.y()));
            bla.replace(1,(particlePosition_.component(0)-ppos.x()));

            //Second order Gaussian
            particleVelocity += pstr/(2.0*pi) * bla/sqr(radius) * (1.0 - exp(-sqr(radius/psize)));
        }

        particlePosition_ = oldParticlePosition + rkCoeffs[iStage]*dt*particleVelocity;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vortexParticleGust::vortexParticleGust(const dictionary& dict,const fvMesh& mesh):
gust(mesh),
dict_(dict),
mesh_(mesh),
R_(readScalar(dict.lookup("radius"))),
strength_(readScalar(dict.lookup("strength"))),
origin_(dict.lookup("origin")),
NvRadial_(readLabel(dict.lookup("NvRadial"))),
NvCirc_(readLabel(dict.lookup("NvCircular"))),
Nv_(NvRadial_*NvCirc_+1),
particlesDistributionRadius_(readScalar(dict.lookup("particlesDistributionRadius"))),
particleStrength_(Nv_,0.0),
particleSize_(Nv_,0.0),
particlePosition_(Nv_,vector::zero),
currentTimeIndex_(-1)
{
    Info << "Gust properties of vortexParticleGust: R = " << R_ <<", Strength = " << strength_ <<", origin = " << origin_ <<  endl;

    //Set initial position and strength
    initializeParticles();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vortexParticleGust::~vortexParticleGust()
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const surfaceVectorField& vortexParticleGust::getGustVelocity(const scalar t,const volVectorField& U)
{
    const scalar& pi = mathematicalConstant::pi;
    Info << "Start getGustVelocity() " << endl;
    //Only move particles once per time step
    if(currentTimeIndex_ < mesh_.time().timeIndex())
    {
        moveParticles(U,mesh_.time().deltaT().value());
        currentTimeIndex_ = mesh_.time().timeIndex();
    }

    Info << "Start field evaluation" << endl;

    //Gather all relevant faces (non empty)
    vectorField faceC(mesh_.nFaces(),vector::zero);
    label counter = 0;
    for(int i=0;i<mesh_.nInternalFaces();i++)
    {
        faceC[counter] = mesh_.faceCentres()[i];
        counter++;
    }

    const surfaceVectorField& faceCentres = mesh_.Cf();
    forAll(particlePosition_,iParticle){
        const scalar& pstr = particleStrength_[iParticle];
        const vector& ppos = particlePosition_[iParticle];
        const scalar& psize = particleSize_[iParticle];
        const vectorField distance = faceCentres.internalField() - ppos;
        //const scalarField radius = mag(faceCentres.internalField() - ppos)+SMALL;

        Ugust_.internalField().component(0) += pstr/(2.0*pi) * (-distance.component(1))/sqr(mag(distance)+SMALL) * (1 - exp(-sqr(mag(distance)+SMALL/psize)));
        Ugust_.internalField().component(1) += pstr/(2.0*pi) * (distance.component(0))/sqr(mag(distance)+SMALL) * (1 - exp(-sqr(mag(distance)+SMALL/psize)));

        //bla.replace(0,-(faceC.component(1)-ppos.y()));
        //bla.replace(1,(faceC.component(0)-ppos.x()));
        //Second order Gaussian
        //Ugust_ += pstr/(2.0*pi) * bla/sqr(radius) * (1 - exp(-sqr(radius/psize)));
    }

    /*forAll(mesh_.boundaryMesh(),iPatch)
    {
        if(mesh_.boundaryMesh().types()[iPatch] != "empty")
        {
            const scalar& startInd = mesh_.boundaryMesh()[iPatch].start();
            const label& patchSize = mesh_.boundaryMesh()[iPatch].size();
            for(int i=0;i<patchSize;i++)
            {
                faceC[counter] = mesh_.faceCentres()[startInd+i];
                counter++;
            }
        }
    }
    faceC.setSize(counter);

    //Calculate gust velocity at all faces in one go
    vectorField UgustTmp(counter,vector::zero);
    vectorField bla(counter,vector::zero);
    forAll(particlePosition_,iParticle)
    {
        const scalar& pstr = particleStrength_[iParticle];
        const vector& ppos = particlePosition_[iParticle];
        const scalar& psize = particleSize_[iParticle];
        const scalarField radius = mag(faceC - ppos)+SMALL;

        bla.replace(0,-(faceC.component(1)-ppos.y()));
        bla.replace(1,(faceC.component(0)-ppos.x()));
        
        //Second order Gaussian
        UgustTmp += pstr/(2.0*pi) * bla/sqr(radius) * (1 - exp(-sqr(radius/psize)));
    }

    //Set internal field
    Ugust_.internalField() = SubList<vector>(UgustTmp,mesh_.nInternalFaces(),0);
    Info << "max(Ugust) = " << max(mag(Ugust_.internalField())) << endl;

    //Do all boundary faces
    label startInd = mesh_.nInternalFaces();
    forAll(Ugust_.boundaryField(),iBoundary){
        if(mesh_.boundaryMesh().types()[iBoundary] != "empty")
        {
            const label& patchSize = mesh_.boundaryMesh()[iBoundary].size();
            Ugust_.boundaryField()[iBoundary] == vectorField(SubList<vector>(UgustTmp,patchSize,startInd));
            startInd += patchSize;
        }
        
    }*/

    return Ugust_;
}

const volVectorField& vortexParticleGust::getGustVelocityCells(const scalar t,const volVectorField& U)
{
    return UgustCell_;
}

const bool vortexParticleGust::isActive() const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
