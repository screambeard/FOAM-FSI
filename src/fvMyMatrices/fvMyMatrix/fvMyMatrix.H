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

Class
    Foam::fvMyMatrix

Description
    Finite-Volume matrix.

SourceFiles
    fvMyMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef fvMyMatrix_H
#define fvMyMatrix_H

#include "volFields.H"
#include "surfaceFields.H"
#include "lduMatrix.H"
#include "tmp.H"
#include "autoPtr.H"
#include "dimensionedTypes.H"
#include "zeroField.H"
#include "className.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

/*---------------------------------------------------------------------------*\
                           Class fvMyMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvMyMatrix
:
    public fvMatrix<Type>
{
    // Private data

        //Diagonal without relaxation
        Field<scalar> D0_;

        //Diagonal without relaxation
        Field<Type> S0_;

public:

    ClassName("fvMyMatrix");


    // Constructors

        //- Construct given a field to solve for
        fvMyMatrix
        (
            GeometricField<Type, fvPatchField, volMesh>&,
            const dimensionSet&
        );

        //- Construct as copy
        fvMyMatrix(const fvMyMatrix<Type>&);

        //- Construct from fvMatrix
        fvMyMatrix(const fvMatrix<Type>& fvm);

        //- Construct as copy of tmp<fvMyMatrix<Type> > deleting argument
#       ifdef ConstructFromTmp
        fvMyMatrix(const tmp<fvMyMatrix<Type> >&);
#       endif

        //Construct from tmp<fvMatrix>
        //fvMyMatrix(const tmp<fvMatrix<Type> >&);

        //- Construct from Istream given field to solve for
        fvMyMatrix(GeometricField<Type, fvPatchField, volMesh>&, Istream&);

        //- Clone
        tmp<fvMyMatrix<Type> > clone() const
        {
            return tmp<fvMyMatrix<Type> >
            (
                new fvMyMatrix<Type>(*this)
            );
        }


    // Destructor

        //virtual ~fvMyMatrix();


    // Member functions

        // Access

            //- Relax matrix (for steady-state solution).
            //  alpha = 1 : diagonally equal
            //  alpha < 1 : diagonally dominant
            //  alpha = 0 : do nothing
            void relax(const scalar alpha);

            //- Relax matrix (for steady-state solution).
            //  alpha is read from controlDict
            void relax();

            void undoRelax(const scalar alpha);

            //- Relax matrix (for steady-state solution).
            //  alpha is read from controlDict
            void undoRelax();

            //- Return the matrix scalar diagonal
            tmp<scalarField> D0() const;

            //- Return the central coefficient
            tmp<volScalarField> A0() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

//#include "fvMyMatrix.C"
#ifdef NoRepository
#   include "fvMyMatrix.C"
#endif

#endif

// ************************************************************************* //
