/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


// mtx.cpp header file

// basic numerical matrix class
//

// 2010/07/12 update
// added Frobenius norm
// 2010/08/18 update
// corrected memory leak bug in destroy(): wrong 'if' statement

#ifndef MTX_H
#define MTX_H

#include "Linalg.h"
#include "Vect.h"

#include <iostream>

class CLin_Matrix
{
public:

	//
	// constructors
	//

	// default constructor
	inline CLin_Matrix() : m_iM(0), m_iN(0), m_iMN(0), m_pV(0), m_pRow(0)
	{
	}

	// construct from another already-existing matrix
	inline CLin_Matrix(const CLin_Matrix &A)
	{
		initialize(A.m_iM, A.m_iN);
		copy(A.m_pV);
	}

	// construct matrix of dim 'M' x 'N' and fill it up with value 'value'
	inline CLin_Matrix(CLin_subscript M, CLin_subscript N, double value = 0.0)
	{
		initialize(M,N);
		set(value);
	}

	// construct matrix of dim 'M' x 'N' and initialize with values from matrix of doubles 'v'
	inline CLin_Matrix(CLin_subscript M, CLin_subscript N, const double* v)
	{
		initialize(M,N);
		copy(v);
	}

	//
	// destructor
	//
	inline ~CLin_Matrix()
	{
		destroy();
	}

	//
	// access
	//

	inline CLin_subscript size() const { return m_iMN; }

    inline double *array() { return m_pV; }

    inline const double *array() const { return m_pV; }
	// returns # of rows or cols according to parameter d (1=rows, 2=cols, other=0)
	inline CLin_subscript dim(CLin_subscript d)
	{
#ifdef LINALG_BOUNDS_CHECK
		CLin_assert( d >= 1);
		CLin_assert( d <= 2);
#endif
		if(d==1) {
			return m_iM;
		}
		else if(d==2) {
			return m_iN;
		}
		else {
			return 0;
		}
	}

	inline CLin_subscript num_rows() const { return m_iM; }

	inline CLin_subscript num_cols() const { return m_iN; }

	//
	// methods
	//

	inline CLin_Matrix& newsize(CLin_subscript M, CLin_subscript N)
	{
		if (num_rows() == M && num_cols() == N)
			return *this;

		destroy();
		initialize(M,N);

		return *this;
	}

	//
	// operators
	//

	inline CLin_Matrix& operator=(const CLin_Matrix &A)
	{
		if (m_pV == A.m_pV)
			return *this;

		// no need to re-alloc
		if (m_iM == A.m_iM  && m_iN == A.m_iN)
			copy(A.m_pV);
		else
		{
			destroy();
			initialize(A.m_iM, A.m_iN);
			copy(A.m_pV);
		}

		return *this;
	}

	inline CLin_Matrix& operator=(const double scalar)
	{
		set(scalar);
		return *this;
	}

	inline operator double**(){ return	m_pRow; }

	inline operator double**() const { return m_pRow; }

	inline double* operator[](CLin_subscript i)
	{
#ifdef LINALG_BOUNDS_CHECK
//		CLin_assert(0<=i);
		CLin_assert(i < m_iM) ;
#endif
		return m_pRow[i];
	}
/*
	inline const double* operator[](CLin_subscript i) const
	{
#ifdef LINALG_BOUNDS_CHECK
		CLin_assert(0<=i);
		CLin_assert(i < m_iM) ;
#endif
		return m_pRow[i];
	}
*/
	inline double operator()(CLin_subscript i, CLin_subscript j)
	{
#ifdef LINALG_BOUNDS_CHECK
//		CLin_assert(i>=0);
		CLin_assert(i < m_iM) ;
//		CLin_assert(j>=0);
		CLin_assert(j < m_iN);
#endif
		return	m_pRow[i][j];
	}


#ifdef LINALG_USE_REGIONS

	typedef Region2D<CLin_Matrix<T> > Region;


	Region operator()(const Index1D &I, const Index1D &J)
	{
		return Region(*this, I,J);
	}


	typedef const_Region2D< CLin_Matrix<T> > const_Region;
	const_Region operator()(const Index1D &I, const Index1D &J) const
	{
		return const_Region(*this, I,J);
	}

#endif

	inline void destroy()
	{
		// do nothing, if no memory has been previously allocated
		if (m_pV != NULL) delete [] (m_pV);
		if (m_pRow != NULL) delete [] (m_pRow);
		m_iM = 0;
		m_iN = 0;
		m_iMN = 0;
		m_pV = NULL;
		m_pRow = NULL;
	}

protected:
	// initialize a new matrix (that is, allocate memory and initialize members)
	inline void initialize(CLin_subscript M, CLin_subscript N)
	{
		m_iMN = M*N;
		m_iM = M;
		m_iN = N;

		// allocate memory
		m_pV = new double[m_iMN];
		m_pRow = new double*[M];

		CLin_assert(m_pV  != NULL);
		CLin_assert(m_pRow	!= NULL);

		// initialize vector of pointers to matrix rows
		double* p = m_pV;
		for (CLin_subscript i=0; i<M; i++)
		{
			m_pRow[i] = p;
			p += N ;
		}
	}

	// copy from an array into the matrix
	inline void copy(const double*  v)
	{
		CLin_subscript i;

		for (i=0; i < m_iM * m_iN; i++)
			m_pV[i] = v[i];
	}

	// fill matrix with value 'val'
	inline void set(const double val)
	{
		CLin_subscript i;

		for (i=0; i < m_iM * m_iN; i++)
			m_pV[i] = val;
	}

	// rows
	CLin_subscript m_iM;
	// cols
	CLin_subscript m_iN;
	// total size
	CLin_subscript m_iMN;
	// the actual matrix
	double* m_pV;
	// array of pointers to matrix rows
	double** m_pRow;
};


//	I/O

std::ostream& operator<<(std::ostream &s, const CLin_Matrix &A);
std::istream& operator>>(std::istream &s, CLin_Matrix &A);

// basic matrix operations

inline CLin_Matrix operator+(const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_assert(M==B.num_rows());
	CLin_assert(N==B.num_cols());

	CLin_Matrix tmp(M,N);
	CLin_subscript i,j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			tmp[i][j] = A[i][j] + B[i][j];

	return tmp;
}

inline CLin_Matrix operator-(const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_assert(M==B.num_rows());
	CLin_assert(N==B.num_cols());

	CLin_Matrix tmp(M,N);
	CLin_subscript i,j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			tmp[i][j] = A[i][j] - B[i][j];

	return tmp;
}

inline CLin_Matrix mult_element(const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_assert(M==B.num_rows());
	CLin_assert(N==B.num_cols());

	CLin_Matrix tmp(M,N);
	CLin_subscript i,j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			tmp[i][j] = A[i][j] * B[i][j];

	return tmp;
}

inline CLin_Matrix transpose(const CLin_Matrix &A)
{
	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_Matrix S(N,M);
	CLin_subscript i, j;

	for (i=0; i<M; i++)
		for (j=0; j<N; j++)
			S[j][i] = A[i][j];

	return S;
}

inline CLin_Vector matmult(const CLin_Matrix  &A, const CLin_Vector &x)
{

#ifdef LINALG_BOUNDS_CHECK
	CLin_assert(A.num_cols() == x.dim());
#endif

	CLin_subscript M, N;

	M = A.num_rows();
	N = A.num_cols();

	CLin_Vector tmp(M);
    
    // Y <- alpha*A*X + beta*Y
    cblas_dgemv(CblasRowMajor,
                CblasNoTrans,
                M,
                N,
                1.0, // Scaling factor for the product of matrix A and vector X
                A.array(),
                M,
                x.array(),
                1,
                0.0, // Scaling factor for vector Y
                tmp.array(),
                1
                );

	return tmp;
}

inline CLin_Matrix matmult(const CLin_Matrix &A, const CLin_Matrix &B)
{

#ifdef LINALG_BOUNDS_CHECK
	CLin_assert(A.num_cols() == B.num_rows());
#endif

	CLin_subscript M, N, K;

	M = A.num_rows();
    N = B.num_cols();
	K = A.num_cols();

	CLin_Matrix tmp(M, N);

    // C <- alpha*A*B + beta*C
    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                M,
                N,
                K,
                1.0,      // ALPHA, scaling factor for the product of matrices A and B.
                A.array(),
                M,
                B.array(),
                N,
                0.0,      // BETA, Scaling factor for matrix C
                tmp.array(),
                M);

	return tmp;
}

inline int matmult(CLin_Matrix &C, const CLin_Matrix &A, const CLin_Matrix &B)
{
	CLin_subscript M, N, K;

	CLin_assert(A.num_cols() == B.num_rows());

	M = A.num_rows();
    N = B.num_cols();
	K = A.num_cols();

	C.newsize(M, N);

    // C <- alpha*A*B + beta*C
    cblas_dgemm(CblasRowMajor,
                CblasNoTrans,
                CblasNoTrans,
                M,
                N,
                K,
                1.0,      // ALPHA, scaling factor for the product of matrices A and B.
                A.array(),
                M,
                B.array(),
                N,
                0.0,      // BETA, Scaling factor for matrix C
                C.array(),
                M);

	return 0;
}

inline CLin_Vector operator*(const CLin_Matrix	&A, const CLin_Vector &x)
{
	return matmult(A,x);
}

inline CLin_Matrix operator*(const CLin_Matrix &A, const CLin_Matrix &B)
{
	return matmult(A,B);
}

extern "C" {
    double dlange_(const char * _Nonnull norm,
                   const __LAPACK_int * _Nonnull m,
                   const __LAPACK_int * _Nonnull n,
                   const double * _Nullable a,
                   const __LAPACK_int * _Nonnull lda,
                   double * _Nullable work);
}

// Frobenius norm
inline double FNorm(const CLin_Matrix &A)
{
	double norm;

	long M = A.num_rows();
    long N = A.num_cols();

    norm = dlange_(
        "f",  // norm type
        &M,
        &N,
        A.array(),
        &M,
        NULL
    );
	return norm;
}

// Frobenius norm, complex matrix
inline double FNorm(const CLin_Matrix &ARe, const CLin_Matrix &AIm)
{
	CLin_subscript M, N;
	CLin_subscript i, j;
	const double *rowiRe, *rowiIm;
	double norm;

	M = ARe.num_rows();
	N = ARe.num_cols();

	CLin_assert(AIm.num_rows() == M);
	CLin_assert(AIm.num_cols() == N);

	if(AIm.num_rows() != M || AIm.num_cols() != N) {
		// impossible value to signal an issue (norm should be positive)
		return -1.0;
	}

	norm = 0;
	for (i=0; i<M; i++) {
		rowiRe = ARe[i];
		rowiIm = AIm[i];
		for (j=0; j<N; j++) {
			norm += rowiRe[j]*rowiRe[j] + rowiIm[j]*rowiIm[j];
		}
	}

	norm = sqrt(norm);

	return norm;
}


#endif // MTX_H
