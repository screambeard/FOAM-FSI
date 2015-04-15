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
    A special matrix type and solver, designed for finite volume
    solutions of scalar equations.
    Face addressing is used to make all matrix assembly
    and solution loops vectorise.

\*---------------------------------------------------------------------------*/

#ifndef fvMyMatricesFwd_H
#define fvMyMatricesFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
class fvMyMatrix;

typedef fvMyMatrix<scalar> fvScalarMyMatrix;
typedef fvMyMatrix<vector> fvVectorMyMatrix;
typedef fvMyMatrix<sphericalTensor> fvSphericalTensorMyMatrix;
typedef fvMyMatrix<symmTensor> fvSymmTensorMyMatrix;
typedef fvMyMatrix<symmTensor4thOrder> fvSymmTensor4thOrderMyMatrix;
typedef fvMyMatrix<diagTensor> fvDiagTensorMyMatrix;
typedef fvMyMatrix<tensor> fvTensorMyMatrix;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
