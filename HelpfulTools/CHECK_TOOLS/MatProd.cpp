// OpenMP optimized matrix mutiplication 
// (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

// to run locally ##########################################################
// g++ -O2 -std=c++11 MatProd.cpp -lm -fopenmp

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

#define NATOM 2884                                                     // Dimesnion of Matrix                                                   // Order of Taylor Expansion

#ifndef NO_OMP
    #include <omp.h>
#endif


typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);


inline int fq(int i, int j, int N)
/*
 gives back MAT[i,j] of vector with leading dimension N
 */
{
	return i*N+j;
}

// Random Complex Number creator
std::complex<double> random_complex()
{
    static std::mt19937 rng( std::time(nullptr) ) ;
    static std::uniform_real_distribution<double> distr( -1.0, 1.0 ) ;

    return { distr(rng), distr(rng) } ;
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
void times1(Vec &A, Vec &B, Vec &C)
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


template <class Vec>
void times2(Vec &A, Vec &B, Vec &C)
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
					C[fq(i,j,NATOM)] += A[fq(i,k,NATOM)]*B[fq(k,j,NATOM)];
				}
			}
}


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


void diagonalize(cvec &Hk, dvec &evals)
{
/*
 Diagonalization of matrix Hk (NATOM x NATOM). Writes eiegenvalues to vector evals and eigenvectors (normalized!) to matrix Hk
 */
	zheev_(&jobz, &uplo, &matsize, &Hk[0], &lda, &evals[0], &work[0], &lwork, &rwork[0], &info);
	assert(!info);
}

int main (){
	//************** OPEN_MP INIT **************************************
#ifndef NO_OMP 	  
	cout << "# of processes " << omp_get_num_procs() << endl;
#pragma omp parallel 
	cout << "Thread " << omp_get_thread_num() << " out of " << omp_get_num_threads() << " says hello!" << endl;     
#endif
	//******************************************************************
	
	double dtime0, dtime1, dtime11, dtime2, dtime3, dtime4;
	cvec *MATRIX3;
	dvec evals(NATOM);
	
	// create random Hermirian matrix	
	cvec MATRIX1(NATOM*NATOM);
	cvec MATRIX2(NATOM*NATOM);
	MATRIX3 = new cvec(NATOM*NATOM);
	cvec MATRIX4(NATOM*NATOM);
	cvec RESULT1(NATOM*NATOM);
  
	for(int i=0 ; i<NATOM; ++i ) 
    {
		for(int j=i ; j < NATOM; ++j )
		{
			const cdouble value = random_complex();
			if(i==j) MATRIX1[fq(i,i,NATOM)] = value.real(); 	// diagonal elements must be real
			else
			{
				MATRIX1[fq(i,j,NATOM)] = value;
				MATRIX1[fq(j,i,NATOM)] = conj(value); 		// complex conjugate of value
			}
		}  
    }
    cout << endl << "Random Hermitian Matrix: " << endl;
    for(int i=0 ; i<NATOM; ++i )
	{	
		for(int j=0 ; j<NATOM; ++j )
		{
			MATRIX2[fq(i,j,NATOM)] =  MATRIX1[fq(i,j,NATOM)];
			MATRIX4[fq(i,j,NATOM)] =  MATRIX1[fq(i,j,NATOM)];
			(*MATRIX3)[fq(i,j,NATOM)] =  MATRIX1[fq(i,j,NATOM)];
		}
	}
	
	//Diagonalization
	dtime0 = omp_get_wtime();
	const clock_t begin_time0 = clock();	
	diagonalize(MATRIX4, evals);
	dtime0 = omp_get_wtime() - dtime0;
    cout << endl << "Diagonalization()" << endl;
    cout << "Calculations time (MPI): " << float(clock() - begin_time0)/CLOCKS_PER_SEC << " seconds" << endl; 
    cout << "Calculations time (OMP): " << dtime0 << " seconds" << endl;
	
	//times1	
	dtime1 = omp_get_wtime();
	const clock_t begin_time1 = clock();
    times1(MATRIX1, MATRIX2, RESULT1);
    dtime1 = omp_get_wtime() - dtime1;
    cout << endl << "times+transpose()" << endl;
    cout << "Calculations time (MPI): " << float(clock() - begin_time1)/CLOCKS_PER_SEC << " seconds" << endl; 
    cout << "Calculations time (OMP): " << dtime1 << " seconds" << endl;
    
    //times1 HEAP	
	dtime11 = omp_get_wtime();
	const clock_t begin_time11 = clock();
    times1(MATRIX1, MATRIX3[0], RESULT1);
    dtime11 = omp_get_wtime() - dtime11;
    cout << endl << "times+transpose() (with pointer)" << endl;
    cout << "Calculations time (MPI): " << float(clock() - begin_time11)/CLOCKS_PER_SEC << " seconds" << endl; 
    cout << "Calculations time (OMP): " << dtime11 << " seconds" << endl;
    
    
    //times2
    dtime2 = omp_get_wtime();
	const clock_t begin_time2 = clock();
    times2(MATRIX1, MATRIX2, RESULT1);
    dtime2 = omp_get_wtime() - dtime2;
    cout << endl << "times WITHOUT transpose()" << endl;
    cout << "Calculations time (MPI): " << float(clock() - begin_time2)/CLOCKS_PER_SEC << " seconds" << endl; 
    cout << "Calculations time (OMP): " << dtime2 << " seconds" << endl;
    
    //times_dn
    dtime3 = omp_get_wtime();
	const clock_t begin_time3 = clock();
    times_dn(MATRIX1, MATRIX2, RESULT1);
    dtime3 = omp_get_wtime() - dtime3;
    cout << endl << "times_dn" << endl;
    cout << "Calculations time (MPI): " << float(clock() - begin_time3)/CLOCKS_PER_SEC << " seconds" << endl; 
    cout << "Calculations time (OMP): " << dtime3 << " seconds" << endl;
    
    //times_nd
    dtime4 = omp_get_wtime();
	const clock_t begin_time4 = clock();
    times_nd(MATRIX1, MATRIX2, RESULT1);
    dtime4 = omp_get_wtime() - dtime4;
    cout << endl << "times_nd" << endl;
    cout << "Calculations time (MPI): " << float(clock() - begin_time4)/CLOCKS_PER_SEC << " seconds" << endl; 
    cout << "Calculations time (OMP): " << dtime4 << " seconds" << endl;

	delete MATRIX3;
}

