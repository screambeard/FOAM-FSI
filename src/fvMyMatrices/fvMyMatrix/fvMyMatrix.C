/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvPatchFields.H"
#include "transformField.H"

//defineTypeNameAndDebug(fvMyMatrix, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fvMyMatrix<Type>::fvMyMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds
)
:
    fvMatrix<Type>(psi,ds),
    D0_(psi.size(), pTraits<scalar>::zero),
    S0_(psi.size(), pTraits<Type>::zero)
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>(GeometricField<Type, fvPatchField, volMesh>&,"
               " const dimensionSet&) : "
               "constructing fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }
}


template<class Type>
Foam::fvMyMatrix<Type>::fvMyMatrix(const fvMyMatrix<Type>& fvm)
:
    fvMatrix<Type>(fvm),
    D0_(fvm.D0_),
    S0_(fvm.S0_)
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>::fvMyMatrix(const fvMyMatrix<Type>&) : "
            << "copying fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }
}


template<class Type>
Foam::fvMyMatrix<Type>::fvMyMatrix(const fvMatrix<Type>& fvm)
:
    fvMatrix<Type>(fvm),
    D0_(fvm.D()),
    S0_(fvm.source())
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>::fvMyMatrix(const fvMatrix<Type>&) : "
            << "copying fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }
}

#ifdef ConstructFromTmp
template<class Type>
Foam::fvMyMatrix<Type>::fvMyMatrix(const tmp<fvMyMatrix<Type> >& tfvm)
:
    fvMatrix<Type>
    (
        const_cast<fvMyMatrix<Type>&>(tfvm()),
        tfvm.isTmp()
    ),
    D0_
    (
        const_cast<fvMyMatrix<Type>&>(tfvm()).D0_,
        tfvm.isTmp()
    ),
    S0_
    (
        const_cast<fvMyMatrix<Type>&>(tfvm()).S0_,
        tfvm.isTmp()
    )
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>::fvMyMatrix(const tmp<fvMyMatrix<Type> >&) : "
            << "copying fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }

    tfvm.clear();
}
#endif

/*template<class Type>
Foam::fvMyMatrix<Type>::fvMyMatrix(const tmp<fvMatrix<Type> >& tfvm)
:
    fvMatrix<Type>
    (
        tfvm()
    ),
    D0_
    (
        this->D()
    ),
    S0_
    (
        this->source()
    )
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>::fvMyMatrix(const tmp<fvMatrix<Type> >&) : "
            << "copying fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }
    Info << "Here I Am" << endl;
}*/

template<class Type>
Foam::fvMyMatrix<Type>::fvMyMatrix
(
    GeometricField<Type, fvPatchField, volMesh>& psi,
    Istream& is
)
:
    fvMatrix<Type>(psi,is),
    D0_(is),
    S0_(is)
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>"
               "(GeometricField<Type, fvPatchField, volMesh>&, Istream&) : "
               "constructing fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }
}


/*template<class Type>
Foam::fvMyMatrix<Type>::~fvMyMatrix()
{
    if (debug)
    {
        Info<< "fvMyMatrix<Type>::~fvMyMatrix<Type>() : "
            << "destroying fvMyMatrix<Type> for field " << this->psi().name()
            << endl;
    }
}*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMyMatrix<Type>::relax(const scalar alpha)
{
    //Not sure how to undo the relaxation. H and source are the same when A and S are set back to previous values.
    //However, A is not the same. Probably to do with all the offdiagonal (boundary) parts it's taking into account.
    //Not a complete class and therefore now commented to prevent misuse.

    //First do relaxation with 1 for storage later
    /*fvMatrix<Type>::relax(1);
    Field<Type>& S = this->source();
    scalarField& D = this->diag();

    //Store the current diagonal and source for later use
    D0_ = D;
    S0_ = S;*/

    fvMatrix<Type>::relax(alpha);
}

template<class Type>
void Foam::fvMyMatrix<Type>::relax()
{
    if (this->psi().mesh().solutionDict().relax(this->psi().name()))
    {
        relax(this->psi().mesh().solutionDict().relaxationFactor(this->psi().name()));
    }
    else
    {
        if (debug)
        {
            InfoIn("void fvMyMatrix<Type>::relax()")
                << "Relaxation factor for field " << this->psi().name()
                << " not found.  Relaxation will not be used." << endl;
        }
    }
}

template<class Type>
void Foam::fvMyMatrix<Type>::undoRelax(const scalar alpha)
{
    //Not sure how to undo the relaxation. H and source are the same when A and S are set back to previous values.
    //However, A is not the same. Probably to do with all the offdiagonal (boundary) parts it's taking into account.
    //Not a complete class and therefore now commented to prevent misuse.
    /*if (alpha <= 0)
    {
        return;
    }

    //Reset diagonal
    scalarField& D = this->diag();
    D = D0_;

    //Reset source
    Field<Type>& S = this->source();
    S = S0_;*/
}


template<class Type>
void Foam::fvMyMatrix<Type>::undoRelax()
{
    if (this->psi().mesh().solutionDict().relax(this->psi().name()))
    {
        undoRelax(this->psi().mesh().solutionDict().relaxationFactor(this->psi().name()));
    }
    else
    {
        if (debug)
        {
            InfoIn("void fvMyMatrix<Type>::undoRelax()")
                << "Relaxation factor for field " << this->psi().name()
                << " not found.  Relaxation will not be used." << endl;
        }
    }
}

template<class Type>
Foam::tmp<Foam::scalarField> Foam::fvMyMatrix<Type>::D0() const
{
    tmp<scalarField> tdiag(new scalarField(D0_));
    this->addCmptAvBoundaryDiag(tdiag());
    return tdiag;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::fvMyMatrix<Type>::A0() const
{
    tmp<volScalarField> tAphi
    (
        new volScalarField
        (
            IOobject
            (
                "A("+this->psi().name()+')',
                this->psi().instance(),
                this->psi().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->psi().mesh(),
            this->dimensions()/this->psi().dimensions()/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    tAphi().internalField() = D0()/this->psi().mesh().V();
    tAphi().correctBoundaryConditions();

    return tAphi;
}

// ************************************************************************* //
