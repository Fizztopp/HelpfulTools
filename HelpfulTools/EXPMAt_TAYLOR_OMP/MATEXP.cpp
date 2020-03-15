// Code to calculate Exp[Matrix] via diagonalization and Taylor Expension: Calculation times and deviation are compared in the end (OpenMP parlallel)
// (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

// to run locally ##########################################################
// g++ -O2 -std=c++11 MATEXP.cpp -llapack -lblas -lm -fopenmp	
// export OMP_NUM_THREADS=4	            
// ./a.out                                         

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
#include <random>
#include <ctime>

using namespace std;

#define NATOM 2000                                                      // Dimesnion of Matrix
#define ORDER 100                                                       // Order of Taylor Expansion

#ifndef NO_OMP
    #include <omp.h>
#endif

//#define TAYLOR_ONLY                                                   // Don't caluylate matrixexponent by diaganoalization (hight numerical effort for high dimension!)

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
int     matsize = NATOM;      										    // The order of the matrix A.  N >= 0
int     lda = NATOM;            										// The leading dimension of the array A.  LDA >= max(1,N)
int     lwork = 2*NATOM-1;      										// The length of the array WORK.  LWORK >= max(1,2*N-1)
double  rwork[3*NATOM-2];       										// dimension (max(1, 3*N-2))
cdouble work[2*NATOM-1];        										// dimension (MAX(1,LWORK)) On exit, if INFO = 0, WORK(1) returns the optimal LWORK
int	    info;

// REMEMBER: Pointer points to an address which can be empty, call by reference cannot.
// Call by reference uses memory of stack, which is limited but acess ist faster. After programm is finished memory is cleared automatically. Memory which is allocated by
// pointers is part of the heap. This is quasi unlimited but slower. After use memory has do be freed by hand -> delete()! (memory leaks)

// Random Complex Number creator
std::complex<double> random_complex()
{
    static std::mt19937 rng( std::time(nullptr) ) ;
    static std::uniform_real_distribution<double> distr( -1.0, 1.0 ) ;

    return { distr(rng), distr(rng) } ;
}


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
void PRINT(Vec &A)
{
	int N = int(abs(sqrt(A.size())));
    for(int i=0 ; i<N; ++i )
	{	
		for(int j=0 ; j<N; ++j )
		{
			cout << A[fq(i,j,N)] << " ";
		}
		cout << endl;
	}	
}


template <class Vec>
inline double NORM(Vec &A)
{
	double norm=0.;
	int N = int(abs(sqrt(A.size())));
    for(int i=0; i<N*N; ++i)
		norm += real(A[i]*conj(A[i]));
	return sqrt(norm);	
}


template <class Vec>
dvec ELEMENTWISE(Vec &A, Vec &B)
{
	double diff;
	dvec diff_MAX = {0.0, 0.0};
	int N = int(abs(sqrt(A.size())));
    for(int i=0; i<N*N; ++i)
    {
		diff = sqrt(real(A[i]*conj(A[i]))) - sqrt(real(B[i]*conj(B[i])));
		if(diff>diff_MAX[0]){
			diff_MAX[0] = diff;
			diff_MAX[1] = double(i);	
		}	
	}
	return diff_MAX;	
}


template <class Vec>
void times(Vec &A, Vec &B, Vec &C)
{
	Vec TEMP(NATOM*NATOM);
	for(int i=0; i<NATOM; i++) {
		for(int j=0; j<NATOM; j++) {
			TEMP[fq(j,i,NATOM)] = B[fq(i,j,NATOM)];
		}
	}
#ifndef NO_OMP 		
		#pragma omp parallel for
#endif 		
		for(int i=0; i<NATOM; ++i)
			for(int j=0; j<NATOM; ++j)
			{
				C[fq(i,j,NATOM)] = 0.;
				for(int k=0; k<NATOM; ++k)
				{
					C[fq(i,j,NATOM)] += A[fq(i,k,NATOM)]*TEMP[fq(j,k,NATOM)];
				}
			}
}


template <class Vec>
void times_dn(Vec &A, Vec &B, Vec &C)
{
	Vec TEMP1(NATOM*NATOM);
	Vec TEMP2(NATOM*NATOM);
	for(int i=0; i<NATOM; i++) {
		for(int j=0; j<NATOM; j++) {
			TEMP1[fq(j,i,NATOM)] = A[fq(i,j,NATOM)];
			TEMP2[fq(j,i,NATOM)] = B[fq(i,j,NATOM)];
		}
	}	
#ifndef NO_OMP 		
		#pragma omp parallel for
#endif 		
		for(int i=0; i<NATOM; ++i)
			for(int j=0; j<NATOM; ++j)
			{
				C[fq(i,j,NATOM)] = 0.;
				for(int k=0; k<NATOM; ++k)
				{
					C[fq(i,j,NATOM)] += conj(TEMP1[fq(i,k,NATOM)])*TEMP2[fq(j,k,NATOM)];
				}
			}        
}


template <class Vec>
void times_nd(Vec &A, Vec &B, Vec &C)
{
#ifndef NO_OMP 		
		#pragma omp parallel for
#endif 		
		for(int i=0; i<NATOM; ++i)
			for(int j=0; j<NATOM; ++j)
			{
				C[fq(i,j,NATOM)] = 0.;
				for(int k=0; k<NATOM; ++k)
				{
					C[fq(i,j,NATOM)] += A[fq(i,k,NATOM)]*conj(B[fq(j,k,NATOM)]);
				}
			}
}


inline double delta(int a, int b)
/*
 Delta function
 */
{
	if (a==b)
		return 1.;
	else
		return 0.;
}


void ExpMat_Taylor(cvec &A, cvec *TEMP1, cvec *TEMP2, cvec &RESULT)
{
	if(ORDER==0)
	{ 
		for(int i=0; i<NATOM; ++i) 
			for(int j=0; j<NATOM; ++j)
				RESULT[fq(i,j,NATOM)] = 1.0*delta(i,j);
	}	
	else if(ORDER==1)
	{
		for(int i=0; i<NATOM; ++i) 
			for(int j=0; j<NATOM; ++j)
				RESULT[fq(i,j,NATOM)] = A[fq(i,j,NATOM)]+delta(i,j); 
	}
	else
	{
		double fac = 1.0;
		for(int i=0; i<NATOM; ++i) 
			for(int j=0; j<NATOM; ++j)
			{
				RESULT[fq(i,j,NATOM)] = A[fq(i,j,NATOM)]+delta(i,j); 
				(*TEMP1)[fq(i,j,NATOM)] = A[fq(i,j,NATOM)];	 
			}
		
		for(int order=2; order<=ORDER; ++order)
		{
			fac *= 1./double(order);
			//cout << "Order: " <<  order << endl;	
			times(TEMP1[0], A, TEMP2[0]);
			TEMP1[0] = TEMP2[0];                                
				//cout << "fac: " <<  fac << endl;
			for(int i=0; i<NATOM*NATOM; ++i) 
				RESULT[i] += fac*(*TEMP1)[i]; 		
		}
	}			
}

int main (){
	//************** OPEN_MP INIT **************************************
#ifndef NO_OMP 	  
	cout << "# of processes " << omp_get_num_procs() << endl;
#pragma omp parallel 
	cout << "Thread " << omp_get_thread_num() << " out of " << omp_get_num_threads() << " says hello!" << endl;     
#endif
	//******************************************************************
	
	double dtime;
	
	// create random Hermirian matrix	
	cvec MATRIX(NATOM*NATOM);
	cvec MATRIX_COPY(NATOM*NATOM);
	cvec *TEMP1 = new cvec(NATOM*NATOM);
	cvec *TEMP2 = new cvec(NATOM*NATOM);
	cvec RESULT1(NATOM*NATOM);
	cvec RESULT2(NATOM*NATOM);
	dvec evals(NATOM); 
  
	for(int i=0 ; i<NATOM; ++i ) 
    {
		for(int j=i ; j < NATOM; ++j )
		{
			const cdouble value = random_complex();
			if(i==j) MATRIX[fq(i,i,NATOM)] = value.real(); 	// diagonal elements must be real
			else
			{
				MATRIX[fq(i,j,NATOM)] = value;
				MATRIX[fq(j,i,NATOM)] = conj(value); 		// complex conjugate of value
			}
		}  
    }
    cout << endl << "Random Hermitian Matrix: " << endl;
    for(int i=0 ; i<NATOM; ++i )
	{	
		for(int j=0 ; j<NATOM; ++j )
		{
			MATRIX_COPY[fq(i,j,NATOM)] =  MATRIX[fq(i,j,NATOM)];
		}
	}
		
	//PRINT(MATRIX);
	
	//Calculate Matrix Exponential by diagonalization
#ifndef	TAYLOR_ONLY
	const clock_t begin_time1 = clock();
	diagonalize(MATRIX_COPY, evals); 
	for(int i=0; i<NATOM; i++)
	{
		(*TEMP1)[fq(i,i,NATOM)] = exp(evals[i]);
	}
	times(TEMP1[0], MATRIX_COPY, TEMP2[0]);                                         // S @ H @ S^(-1) = H_D --> nk = S^(-1) @ nk_D @ S
	times_dn(MATRIX_COPY, TEMP2[0], RESULT1);	
	
	cout << endl << "Matrix Exponential by diaglonalization: " << endl;
	cout << "Calculations time: " << float(clock() - begin_time1)/CLOCKS_PER_SEC << " seconds" << endl; 
	//PRINT(RESULT1);
#endif
	
	dtime = omp_get_wtime();
	//Calculate Matrix Exponential by Taylor Series
	ExpMat_Taylor(MATRIX, TEMP1, TEMP2, RESULT2);
	dtime = omp_get_wtime() - dtime;
	cout << endl << "Matrix Exponential by Taylor Series: " << endl;
	cout << "Calculations time (OMP): " << dtime << " seconds" << endl; 
	//PRINT(RESULT2);
	
		
	cout << endl << "Dimension of matrix: " << NATOM << "x" << NATOM << endl;
	cout << "Expansion Order: " << ORDER << endl;
#ifndef	TAYLOR_ONLY	
	cout << "Abs. deviation of norm: " << NORM(RESULT1)-NORM(RESULT2) << endl;
	int THREAD = omp_get_thread_num();
	dvec MAX_DIFF = ELEMENTWISE(RESULT1, RESULT2);
	if(THREAD==0) cout << "max elementwise deviation: " << MAX_DIFF[0] << " at position: " << MAX_DIFF[1] << endl;
#endif	

	delete TEMP1, TEMP2;
}

