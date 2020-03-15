//run locally ##########################################################
// g++ -O2 -std=c++11 LAPACK_EXAMPLE.cpp -llapack -lblas -lm		    -> compile with c++ (#define NO_MPI)
// ./a.out
// mpic++ -O2 -std=c++11 LAPACK_EXAMPL.cpp -llapack -lblas -lm      -> compile with MPI
// mpiexec -n 4 a.out                                                   -> run parallel on 4 cores

//run on DRACo ###########################################################

//export LD_LIBRARY_PATH="$MKL_HOME/lib/intel64"
//module load mkl gcc impi
//mpicxx -O2 -std=c++11  LAPACK_EXAMPL.cpp -L$MKL_HOME/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm


#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <assert.h>
#include <iterator>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

#define NORB 1


using namespace std;

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);

//LAPACK (Fortran 90) functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//routine to find eigensystem of Hk
extern "C" {
// computes the eigenvalues and, optionally, the left and/or right eigenvectors for HE matrices
void zheev_(char* jobz, char* uplo, int* N, cdouble* H, int* LDA, double* W, cdouble* work, int* lwork, double* rwork, int *info);
}

char    jobz = 'V';             									    //'N','V':  Compute eigenvalues only/+eigenvectors
char    uplo = 'U';              										//'U','L':  Upper/Lower triangle of H is stored
int     matsize = 2*NORB;      										    // The order of the matrix A.  N >= 0
int     lda = 2*NORB;            										// The leading dimension of the array A.  LDA >= max(1,N)
int     IPIV[2*NORB];            										// IPIV is INTEGER array, dimension (N)
int     lwork = 2*2*NORB-1;      										// The length of the array WORK.  LWORK >= max(1,2*N-1)
double  rwork[3*2*NORB-2];       										// dimension (max(1, 3*N-2))
cdouble work[2*2*NORB-1];        										// dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK
int	    info;

// REMEMBER: Pointer points to an address which can be empty, call by reference cannot.
// Call by reference uses memory of stack, which is limited but acess ist faster. After programm is finished memory is cleared automatically. Memory which is allocated by
// pointers is part of the heap. This is quasi unlimited but slower. After use memory has do be freed by hand -> delete()! (memory leaks)


void diagonalize(cvec &Hk, dvec &evals)
{
/*
 Diagonalization of matrix Hk[64]. Writes eiegenvalues to vector evals and eigenvectors to Hk
 */
	zheev_(&jobz, &uplo, &matsize, &Hk[0], &lda, &evals[0], &work[0], &lwork, &rwork[0], &info);
	assert(!info);
}


inline int fq(int i, int j, int N)
/*
 gives back MAT[i,j] of vector with leading dimension N
 */

{
	return i*N+j;
}


template <class Vec>
void transpose(Vec &mat)
/*
 transposes of mat
 */
{ 
	cdouble temp;
	int N = int(abs(sqrt(mat.size())));
	for(int i = 1; i < N; i++)
	{
		for (int j = 0; j < i; j++) 
		{
         temp = mat[fq(i,j,N)];
         mat[fq(i,j,N)] = mat[fq(j,i,N)];
         mat[fq(j,i,N)] = temp;
        }
	}
}


void conj_T(cvec &mat)
/*
 conjugate-transposes of mat
 */
{ 
	cdouble temp;
	int N = int(abs(sqrt(mat.size())));
	for(int i = 1; i < N; i++)
	{
		for (int j = 0; j < i; j++) 
		{
         temp = mat[fq(i,j,N)];
         mat[fq(i,j,N)] = mat[fq(j,i,N)];
         mat[fq(j,i,N)] = temp;
        }
	}
	for(int i=0; i<N*N; i++)
		mat[i] = conj(mat[i]);
}


template <class Vec>
void times(Vec &A, Vec &B, Vec &C)
{
	int N = int(abs(sqrt(A.size())));
    for(int i=0; i<N; ++i)
        for(int j=0; j<N; ++j)
        {
            C[fq(i,j,N)] = 0.;
            for(int k=0; k<N; ++k)
            {
                C[fq(i,j,N)] += A[fq(i,k,N)]*B[fq(k,j,N)];
            }
        }
}


template <class Vec>
void times_dn(Vec &A, Vec &B, Vec &C)
{
	int N = int(abs(sqrt(A.size())));
	for(int i=0; i<N; ++i)
        for(int j=0; j<N; ++j)
        {
            C[fq(i,j,N)] = 0.;
            for(int k=0; k<N; ++k)
            {
                C[fq(i,j,N)] += conj(A[fq(k,i,N)])*B[fq(k,j,N)];
            }
        }
}


template <class Vec>
void times_nd(Vec &A, Vec &B, Vec &C)
{
	int N = int(abs(sqrt(A.size())));
	for(int i=0; i<N; ++i)
        for(int j=0; j<N; ++j)
        {
            C[fq(i,j,N)] = 0.;
            for(int k=0; k<N; ++k)
            {
                C[fq(i,j,N)] += A[fq(i,k,N)]*conj(B[fq(j,k,N)]);
            }
        }
}


template <class Vec>
void mat_prod(Vec &mat_return, Vec &mat1, Vec &mat2)
/*
 Matrix product of quadratic matrices: mat_return = mat1 @ mat2:
 */
{
	int N = int(abs(sqrt(mat1.size())));
    for(int i=0; i<N; i++)
    {
    	for(int j=0; j<N; j++)
    	{
    		mat_return[fq(i,j,N)] = 0.;
    		for(int k=0; k<N; k++)
    		{
    			mat_return[fq(i,j,N)] += mat1[fq(i,k,N)] * mat2[fq(k,j,N)];
    		}
    	}
    }
}


int main (){
	
	cvec temp(4,0.);
	cvec SI(4,0.);
	cvec Hk(4); //= {1.,2.+II,2.-II,3.};
	Hk[fq(0,0,2)] = 1.;
	Hk[fq(0,1,2)] = 2.+II;
	Hk[fq(1,0,2)] = 2.-II;
	Hk[fq(1,1,2)] = 3.;
	
	dvec evals(NORB, 0.);
    
    cout << "Hamiltonian --------------------------" << endl;
	for (int i=0; i<4; i++)
	{
		cout<<Hk[i]<<endl;
	}
	
	diagonalize(Hk, evals);
	
	cout << endl;
	cout << "EigenValues --------------------------" << endl;
	for (int i=0; i<2; i++)
	{
		cout<< " " << evals[i] << endl;
	}
	
	SI = Hk;
	
	Hk[fq(0,0,2)] = 1.;
	Hk[fq(0,1,2)] = 2.+II;
	Hk[fq(1,0,2)] = 2.-II;
	Hk[fq(1,1,2)] = 3.;
	
	times_nd(Hk, SI, temp);
	times(SI, temp, Hk);												// S @ H @ S^(-1) = H_D 
	
	
	cout << endl;														            
	cout << "diagonalized Hamiltonian --------------" << endl;
	
	for (int i=0; i<4; i++)
	{
		cout<<Hk[i]<<endl;
	}
	
	cout << endl;
	cout << "original Hamiltonian ------------------" << endl;
	
	times(Hk, SI, temp);                                                // S @ H @ S^(-1) = H_D --> nk = S^(-1) @ nk_D @ S
	times_dn(SI, temp, Hk);
	
	for (int i=0; i<4; i++)
	{
		cout<<Hk[i]<<endl;
	}

	
}

