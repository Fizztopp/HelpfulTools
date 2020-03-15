// This codes calculates from a given List of vectors (here from a 3d BZ of Pyrochlore ["k_BZ_full.txt"]) all possible scattering quadruplets which fullfill momentum conservation. 
// For MPI parallelization two loops are necessary. 1st loop calcultes number of quadruplets per process, 2nd loop calculates actually calculates quadruplets.  
// (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

//run locally without MPI ##########################################################

// MPI ##############################################################
// g++ -O2 -std=c++11 BOLTZMANN_QUADRUPLETS.cc -lm -DNO_MPI 
// ./a.out 

//run locally with MPI ##########################################################

// MPI ##############################################################
// mpic++ -O2 -std=c++11 BOLTZMANN_QUADRUPLETS.cc -llapack -lblas -lm  
// mpiexec -n 4 a.out  

//run on DRACO ##########################################################

//#!/bin/bash
//# Standard output and error:
//#SBATCH -o ./examplejob.out
//#SBATCH -e ./examplejob.err
//#SBATCH -D ./
//#SBATCH --nodes=8
//#SBATCH --ntasks-per-node=32
//#SBATCH --cpus-per-task=1
//## Request 500 GB of main Memory per node in Units of MB:
//##SBATCH --mem=512000
//#SBATCH -J ZERO_DRIVE
//#SBATCH --mail-type=none
//#SBATCH --partition=express
//#SBATCH --time=00:30:00

//module load impi
//module load mkl

//export LD_LIBRARY_PATH="$MKL_HOME/lib/intel64"

//mpicxx -O2 -std=c++11  BOLTZMANN_QUADRUPLETS.cc -L$MKL_HOME/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm

//srun ./a.out >log




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

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<int> Ivec;
 
cdouble II(0,1);

#ifndef NO_MPI
    #include <mpi.h>
#endif


// PARAMETERS ##########################################################


inline int fq(int i, int j, int N)
/*
 gives back MAT[i,j] of vector with leading dimension N
 */
{
	return i*N+j;
}


void ReadIn(vector<dvec> &MAT, const string& filename)
{
	ifstream in(filename);
	string record;
	if(in.fail()){
		cout << "file" << filename << "could not be found!" << endl;
	}
	while (getline(in, record))
	{
		istringstream is( record );
		dvec row((istream_iterator<double>(is)),
		istream_iterator<double>());
		MAT.push_back(row);
	}
	in.close();
}

// main() function #####################################################

int main(int argc, char *argv[])
{
	//************** MPI INIT ***************************
	const int root = 0;
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
	if(myrank==root) cout << "\n\tProgram running on " << numprocs << " processors." << endl;
	
	//************** OPEN_MP INIT ***************************
	
	//vector of BZ vectors
	vector<dvec> BZ_IRR;
	ReadIn(BZ_IRR, "k_BZ_full.txt");
	if(myrank==root) cout << "irreducible BZ --> " << BZ_IRR.size() << " points" << endl;
	int num_kpoints_BZ = BZ_IRR.size();                                                          // # of k-points in BZ
    int dim_k = BZ_IRR[0].size();																 // dimension of vectors (here 3)
	
	// CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	const clock_t begin_time = clock();	
	
	
	int count_tot = 0;                                                  // total number of quadruplets
	int count_proc = 0;													// number of quadruplets per process
	double treshold = 1e-10;                                            // accuracy for momentum conservation	
	
	
	dvec temp(3);
    
    // 1st loop: calculate # of quadruplets per process
    if(myrank==root) cout << "Start 1st Loop" << endl;
    	
	for(int k1=myrank; k1<num_kpoints_BZ; k1+=numprocs)
	{
		if(myrank==root) cout << "1st Run: "<< k1 << endl;
		for(int k2=0; k2<num_kpoints_BZ; k2++)
		{
			if(k1==k2) 
				continue;
			else
				for(int k3=0; k3<num_kpoints_BZ; k3++)
				{
					if(k2==k3 || k1==k3) 
						continue;
					else	
						for(int k4=0; k4<num_kpoints_BZ; k4++)
						{
							if(k1==k4 || k2==k4 || k3==k4) 
								continue;
							else
								for(int i=0; i<dim_k; i++)
									temp[i] = abs(BZ_IRR[k1][i] + BZ_IRR[k2][i] - BZ_IRR[k3][i] - BZ_IRR[k4][i]);
								if(temp[0]>treshold || temp[1]>treshold || temp[2]>treshold)
									continue;
								count_proc += 1;
						}	
				}
		}
	}
		
	//}
	// Calculate total number of quadruplets and send to all processes
#ifndef NO_MPI		
	MPI_Allreduce(&count_proc, &count_tot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif	
	
	// 2nd loop: every process calculates a number of count_proc quaduplets                                  	
	int count_help = 0;	
	
	vector<int> QUAD_PROC(count_proc*4);
	vector<int> help;
    
    if(myrank==root) cout << "Start 2nd Loop" << endl;
    	
	for(int k1=myrank; k1<num_kpoints_BZ; k1+=numprocs)
	{
		if(myrank==root) cout << "2nd Run: "<< k1 << endl;
		
		for(int k2=0; k2<num_kpoints_BZ; k2++)
		{
			if(k1==k2) 
				continue;
			else
				for(int k3=0; k3<num_kpoints_BZ; k3++)
				{
					if(k2==k3 || k1==k3) 
						continue;
					else	
						for(int k4=0; k4<num_kpoints_BZ; k4++)
						{
							if(k1==k4 || k2==k4 || k3==k4) 
								continue;
							else
								for(int i=0; i<dim_k; i++)
									temp[i] = abs(BZ_IRR[k1][i] + BZ_IRR[k2][i] - BZ_IRR[k3][i] - BZ_IRR[k4][i]); // momentum has to be conserved for all spatial (linear independent) dimensions!
								if(temp[0]>treshold || temp[1]>treshold || temp[2]>treshold)
									continue;
								help = {k1,k2,k3,k4};
								for(int k=0; k<4; k++)
									QUAD_PROC[fq(count_help,k,4)] = help[k];
								count_help += 1;
						}	
				}
		}
	}
	
	int *COUNTS, *displs, *QUAD  = NULL;
	
	if(myrank==root) COUNTS = (int *)malloc(numprocs*sizeof(int));
	
	// Collect counts of differen processes (count2) in COUNTS
#ifndef NO_MPI		
	MPI_Gather(&count_proc, 1, MPI_INT, COUNTS, 1, MPI_INT, root, MPI_COMM_WORLD);	
#endif					
	
	// Check COUNTS
	int test_count = 0;
	if(myrank==root)
	{
		for(int i=0; i<numprocs; i++)
		{
			cout << COUNTS[i] << " ";
			test_count += COUNTS[i]; 
		}
		cout << "--> Sum_count: " << test_count << " " << "CHECK total_count: " << count_tot <<endl;
	}	
	
	// GATHER arrays (quad) together in QUAD 
#ifndef NO_MPI	
	int sendint = count_proc*4;
	
	if(myrank==root)  cout << "sendint: " << sendint <<  " CHECK: " << COUNTS[myrank]*4 << endl; // COUNTS onlly known by root!
	
	if(myrank==root)
	{
		displs = (int *)malloc(numprocs*sizeof(int));
		displs[0]=0; 
		for(int i=1; i<numprocs; i++)
			displs[i]=(displs[i-1]+COUNTS[i-1]*4);
		for(int i=0; i<numprocs; i++)
		{
			COUNTS[i] = 4*COUNTS[i];
			cout << displs[i] << " ";
		}
		QUAD = (int *)malloc(4*count_tot*sizeof(int));		
	}	
	// root gets all calculted quadruplets from all processes	
	MPI_Gatherv(&QUAD_PROC[0], sendint, MPI_INT, QUAD, COUNTS, displs, MPI_INT, root, MPI_COMM_WORLD);
#endif
	
	if(myrank==root) cout << "Calculations lasted: " << float(clock() - begin_time)/CLOCKS_PER_SEC << " seconds" << endl;
	
	// root saves all quadruplets to file
#ifndef NO_MPI		
	if(myrank==root) 
	{
		ofstream myfile ("QUADRUPLETS.txt");
		if (myfile.is_open())
		{
			for(int cc=0; cc<count_tot; cc++)
			{
				for(int kk=0; kk<4; kk++)
					myfile << QUAD[fq(cc,kk,4)]<< " ";
				myfile  << endl;	
			}	
			myfile.close();
		}
		else cout << "Unable to open file" << endl;	
	}	
#else	
	ofstream myfile ("QUADRUPLETS.txt");
	if (myfile.is_open())
	{
		for(int cc=0; cc<count_proc; cc++)
		{
			for(int kk=0; kk<4; kk++)
				myfile << QUAD_PROC[fq(cc,kk,4)]<< " ";
			myfile  << endl;	
		}	
		myfile.close();
	}
	else cout << "Unable to open file" << endl;		
#endif		
	
#ifndef NO_MPI
	MPI_Finalize();
#endif	
	
	
}
