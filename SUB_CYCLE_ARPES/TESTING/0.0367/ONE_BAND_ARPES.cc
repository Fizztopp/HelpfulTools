// TIGHT-BINDING MODEL FOR InSI
// (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

// This code provides a numerical solution of the Floquet driven 1d InSI chain

// run locally #########################################################
// g++ -O2 -std=c++11 ONE_BAND__ARPES.cc -llapack -lblas -lm		       			-> compile with c++ (#define NO_MPI)
// ./a.out
// mpic++ -O2 -std=c++11 ONE_BAND_ARPES.cc -llapack -lblas -lm             		-> compile with MPI
// mpiexec -n 4 a.out                                                   -> run parallel on 4 cores

// run on DRACO ########################################################

// #!/bin/bash
// # Standard output and error:
// #SBATCH -o ./examplejob.out
// #SBATCH -e ./examplejob.err
// #SBATCH -D ./
// #SBATCH --nodes=1
// #SBATCH -J InSi
// #SBATCH --mail-type=none
// #SBATCH --partition=general
// #SBATCH --time=23:59:59

// module load mkl
// module load impi
// module load mkl
// export LD_LIBRARY_PATH="$MKL_HOME/lib/intel64"
// mpicxx -O2 -std=c++11 ONE_BAND__ARPES.cc -L$MKL_HOME/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

// srun -n 32 ./a.out >log


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


// PARAMETERS ##########################################################

#define PI        3.14159265359

// intrinsic parameters
// electronic
#define NN        128                                                   // # of k-points
#define t0        0.2
#define mu        0.4
#define BETA      100.0                     					        // inverse temperature  300 K -> 40, 40 K --> 250

// Peierls driving
#define Ax_peierls    0.0364                               				// Amplitude of Applied Field in x-direction
#define DELAY		 0.0          //- > 0 fs    

// Tr-ARPES
// probe pulse
#define starttime 0.0 
#define endtime   760.0 //- > 500fs
#define TIMESTEPS 1e3                                                 	// # of timesteps used for G<(t,t') -> timesteps/TIMESTEPS should be int!
#define T_PROBE   380.0   // -> 250fs                      				// delay
#define SIGMA_PROBE 65.0  //- > 100 fs @FWHM                            // RMS width
#define OMEGA_PROBE_MIN -0.9                                            // MINIMUM probe frequency
#define OMEGA_PROBE_MAX  +0.1                                           // MAXIMUM probe frequency
#define N_OMEGA_PROBE 128                                               // # of probed frequencies
#define K_PROBE_MIN 0                                                   // MINIMUM probe moementum 
#define K_PROBE_MAX 128                                                 // MINIMUM probe moementum
#define weightcutoff 1e-6 											    // cutoff for Gaussian weight defines window where ARPES signal is collected

#ifndef NO_MPI
    #include <mpi.h>
#endif

using namespace std;

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);

// DEFINITION OF FUNCTIONS #############################################


//INLINE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inline double HEAVISIDE(int t1, int t2)
/*
 THETA function
 */
{
	if (t1 > t2)
		return 1.;
	else
		return 0.;
}


inline int fq(int i, int j, int N)
/*
 gives back MAT[i,j] of vector with leading dimension N
 */

{
	return i*N+j;
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


template <class Vec>
inline void print(Vec vec)
/*
 print out vector
 */
{
	for(int i=0; i<vec.size(); i++)
		{
	    	cout << vec[i] << " ";
	    }	
	cout << endl;
}


inline double Ax_t(double time)
{
/*
 Peierls field for electrons:
 -energy: energy eigenvalue
 -mu: chemical potential
 */
	if(time>DELAY)
		return -Ax_peierls*(time-DELAY);
	else
		return 0.0;
}


inline double fermi(double energy)
{
/*
 Fermi distribution:
 -energy: energy eigenvalue
 -mu: chemical potential
 */
    return 1./(exp((energy)*BETA) + 1.);
}


inline double gauss(double time, double delay, double sigma)
/*
 Gauss distribution function NOT NORMALIEz
 -time: time coordinate
 -delay: mean expectation value
 -sigma: standard deviation 
 */
{
	return 1./(sigma*sqrt(2.*PI))*exp(-0.5*pow((time-delay)/sigma,2.));
}


// VOID FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
/*
 matrix multiplication
 */
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
/*
 C = conj(A).transpose() \cdot B
 */
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
/*
 C = A \cdot conj(B).transpose()
 */
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


double set_Hk(double k, double time)
/*
 Set electronic Hamiltonian
  -k: quasi momentum 
  -u: nuclear displacement
  -Hk: matrix Hamiltonian Hk[k]
  -time: time coordinate
 */
{
	return -2.*t0*cos(k-Ax_t(time)) - mu;
}


double set_Hk0(double k)
/*
 Set electronic Hamiltonian
  -k: quasi momentum 
  -u: nuclear displacement
  -Hk: matrix Hamiltonian Hk[k]
  -time: time coordinate
 */
{
	return -2.*t0*cos(k) - mu;
}

	
void Hk_bands(dvec &BZ, double time, const string& filename)
/*
 Calculate bands of Hk(k) for high symmetry path and store them in "bands.txt":
  -BANDS: vector of vectors containing eigenvalues
  -DENS: array which contains pseudo spins [dim,orbital] -> magnetization
  -Hk: Matrix of dimension [64]
  -evals: vector to store eigenvalues
  -K_PATH: vector of high-symmetry path vectors
  -MAT_BASIS: Basis vectors in array of dimension [4][3]
 */
{
	ofstream myfile (filename);
	if (myfile.is_open())
	{
		for(int k=0; k<NN; k++)
		{
			myfile << set_Hk(BZ[k], time) << " " ;
		}
	myfile.close();
	}
    else cout << "Unable to open file" << endl;
}


double Iphoto_RAMP(int k, dvec &A_INTEGRAL, double omega, dvec &BZ) 
/*
 Calculate Photo current from Tr_Gless()
 */
{
	double h = (endtime-starttime)/TIMESTEPS;
	double Iph;
	dvec TEMP(TIMESTEPS);
	int k_back_shift1;
	int k_back_shift2;
	double shift1;
	double shift2;
	cdouble G_Tr1;
	cdouble G_Tr2;
	
	for(int t2=0; t2<TIMESTEPS; t2++)
	{
		TEMP[t2] = 0.0;
		for(int t1=0; t1<TIMESTEPS-1; t1++)
		{
			if(gauss(double(t1)*h, T_PROBE, SIGMA_PROBE)*gauss(double(t2)*h, T_PROBE, SIGMA_PROBE)*(2.*PI*pow(SIGMA_PROBE, 2))>weightcutoff)
			{
				shift1 = -Ax_peierls*h*double(t1+t2)/2.;
				k_back_shift1 = k+int(round(shift1/(2.*PI)*double(NN)));
				while (k_back_shift1<0)
					k_back_shift1 += NN;
				while(k_back_shift1>NN-1)
					k_back_shift1 -= NN;			
				shift2 = -Ax_peierls*h*double(t1+1+t2)/2.;
				k_back_shift2 = k+int(round(shift2/(2.*PI)*double(NN)));
				while (k_back_shift2<0)
					k_back_shift2 += NN;
				while (k_back_shift2>NN-1)
					k_back_shift2 -= NN;
				G_Tr1 = II*fermi(set_Hk0(BZ[k_back_shift1]))*exp(-2.*II*set_Hk0(BZ[k])/Ax_peierls*sin(Ax_peierls*h*double(t1-t2)/2.));
				G_Tr2 = II*fermi(set_Hk0(BZ[k_back_shift2]))*exp(-2.*II*set_Hk0(BZ[k])/Ax_peierls*sin(Ax_peierls*h*double(t1+1-t2)/2.));
				TEMP[t2] += 0.5*h*imag( gauss(double(t1)*h, T_PROBE, SIGMA_PROBE)*gauss(double(t2)*h, T_PROBE, SIGMA_PROBE)*exp(II*omega*double(t1-t2)*h)*G_Tr1
									   +gauss(double(t1+1)*h, T_PROBE, SIGMA_PROBE)*gauss(double(t2)*h, T_PROBE, SIGMA_PROBE)*exp(II*omega*double(t1+1-t2)*h)*G_Tr2);
			}
		}	  
	}	
	Iph = 0.0;	
	for(int t2=0; t2<TIMESTEPS-1; t2++)
	{
		Iph += 0.5*h*(TEMP[t2]+TEMP[t2+1]);
	}			
	return Iph;
}


void EDC(int k, dvec &A_INTEGRAL, dvec &IPHOTO, int &myrank,  dvec &BZ)
/*
 Calculate of Energy Distribution Curve
 */
{
	dvec OMEGA(N_OMEGA_PROBE);
	for(int w=0; w<N_OMEGA_PROBE; w++)
	{
		OMEGA[w] =  OMEGA_PROBE_MIN + (OMEGA_PROBE_MAX - OMEGA_PROBE_MIN)/N_OMEGA_PROBE*double(w);
	    IPHOTO[w] = Iphoto_RAMP(k, A_INTEGRAL, OMEGA[w], BZ); 
		if(myrank==0) cout << "w#" << w << " = " << OMEGA[w] << endl; 
	}
	if(myrank==0)
	{
		ofstream myfile ("IPHOTO.txt");
		if (myfile.is_open())
		{
			for(int w=0; w<N_OMEGA_PROBE; w++)
			{
				myfile << IPHOTO[w]<<  endl;
			}
			myfile.close();
		}
    else cout << "Unable to open file" << endl;
	}	
}	


void TrARPES(dvec &BZ, dvec &IPHOTO, dvec &ARPES, dvec &A_INTEGRAL, int &numprocs, int &myrank)
/*
 Calculate of energy(w) and angel(k) resolved Photocurrent signal
 */
{
	for(int k=K_PROBE_MIN+myrank; k<K_PROBE_MAX; k+=numprocs)
	{
		EDC(k, A_INTEGRAL, IPHOTO, myrank, BZ);	
		for(int w=0; w<N_OMEGA_PROBE; w++)
		{		
			ARPES[fq(w, k-K_PROBE_MIN, K_PROBE_MAX-K_PROBE_MIN)] = IPHOTO[w];
		}
	}
	MPI_Allreduce(MPI_IN_PLACE, &ARPES[0], N_OMEGA_PROBE*(K_PROBE_MAX-K_PROBE_MIN), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	// Print to file	
	if(myrank==0)
	{
		print(ARPES);
		ofstream myfile ("DATA/ARPES.txt");
		if (myfile.is_open())
		{
			for(int w=0; w<N_OMEGA_PROBE; w++)
			{
				for(int k=0; k<K_PROBE_MAX-K_PROBE_MIN; k++)
				{
					myfile << ARPES[fq(N_OMEGA_PROBE-w-1,k,K_PROBE_MAX-K_PROBE_MIN)] << " " ;
				}
			myfile  << endl;
			}	
			myfile.close();
		}
    else cout << "Unable to open file" << endl;
	}
}
// main() function #####################################################

int main(int argc, char * argv[])
{
    //************** MPI INIT ***************************
  	int numprocs=1, myrank=0, namelen;
    
#ifndef NO_MPI
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	MPI_Init(&argc, &argv);
  	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  	MPI_Get_processor_name(processor_name, &namelen);
    
	cout << "Process " << myrank << " on " << processor_name << " out of " << numprocs << " says hello." << endl;
	MPI_Barrier(MPI_COMM_WORLD);
    
#endif
	if(myrank==0) cout << "\n\tProgram running on " << numprocs << " processors." << endl;
   
	// DECLARATION AND INTITALIZATIO
	
	//vector of high-symmetry path vectors
	dvec BZ(NN);
	for(int ii=0; ii<NN; ii++)
	{
		BZ[ii] = -PI + 2.*PI*double(ii)/double(NN);                        
	}
	
	
	// ARRAYS for ARPES calculations   
	                
	dvec IPHOTO(N_OMEGA_PROBE);
	dvec ARPES(N_OMEGA_PROBE*(K_PROBE_MAX-K_PROBE_MIN));
	
	
	// CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	double h = (endtime-starttime)/TIMESTEPS;
	const clock_t begin_time = clock();	
	
	// calculate driving field
	double Ax;
	if(myrank==0)
	{
		ofstream myfile ("DATA/DRIVING_t.txt");
		if (myfile.is_open())
		{
			for(int t=0; t<TIMESTEPS; t++)
			{			
				Ax = Ax_t(h*double(t));
				myfile << Ax << " "  << endl;	
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	
	// Precalculate Integral values of A(t)
	int timesteps = 1e5;
	double dt = (endtime-starttime)/timesteps;
	dvec A_INTEGRAL(timesteps);
	A_INTEGRAL[0] = 0.0;
	for(int t=0; t<timesteps-1; t++)
	{
		A_INTEGRAL[t+1] = A_INTEGRAL[t] + dt*0.5*(Ax_t(double(t)*dt)+Ax_t(double(t+1)*dt));
		//cout << A_INTEGRAL[t] << endl;
	}	
	if(myrank==0)
	{
		ofstream myfile ("A_INT.txt");
		if (myfile.is_open())
		{
			for(int t=0; t<TIMESTEPS; t++)
			{			
				myfile << A_INTEGRAL[t*timesteps/TIMESTEPS] << " " << A_INTEGRAL[t*timesteps/TIMESTEPS]/(h*double(t)) << endl;	
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}
	
	// caclulate initial band structure
	if(myrank==0) 
	{
		//caclulation of initial bands
		Hk_bands(BZ, 0.0, "bands0.txt");
	}
	
	// Popagate and get toal energy
	//cout << Ax_peierls << "#######################################################################################" << endl; 
	
	//Calculate trARPES
	//if(myrank==0)Tr_Gless(64, BZ, G_Tr, A_INTEGRAL, myrank);
	//if(myrank==0)EDC(64, A_INTEGRAL, G_Tr, IPHOTO, myrank);
	
	TrARPES(BZ, IPHOTO, ARPES, A_INTEGRAL, numprocs, myrank);
	
	//for(int k=K_PROBE_MIN+myrank; k<K_PROBE_MAX; k+=numprocs)
	//{
	//	if(myrank==0) cout << "k = " << k << endl; 
	//	Tr_Gless(k, BZ, G_Tr, myrank);
	//}	
	//MPI_Allreduce(MPI_IN_PLACE, &G_Tr[0], TIMESTEPS*TIMESTEPS*NN, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
	
	//if(myrank==0) Iphoto(64, A_INTEGRAL, -0.6, G_Tr); 
	
	
	if(myrank==0) cout << "Calculations lasted: " << float(clock() - begin_time)/CLOCKS_PER_SEC << " seconds" << endl;
	
#ifndef NO_MPI
	MPI_Finalize();
#endif		
}

