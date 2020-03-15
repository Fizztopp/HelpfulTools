// mpic++ -O2 MPI_EXAMPLE.cc                							    -> compile with MPI
// mpiexec -n 4 a.out                                                   -> run parallel with 4 processes

#include <vector>
#include <math.h>
#include <algorithm>

//#define NO_MPI                                                        <- switch MPI on/off        										

#ifndef NO_MPI
    #include <mpi.h>
#endif

using namespace std;

typedef vector<double> dvec;                     					    // vectors with real double values



// main() function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main(int argc, char * argv[])
{
    //************** MPI HEADER ***************************
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


//	%%%%%%%%%%%%%%%% ALLOCATION 
	
	int size = 100000;
	dvec vec_test(size);          
																		

// %%%%%%%%%%%%%%%%% CACLULATIONS

	for(int i=myrank; i<size; i+=numprocs)								// each processes will write the value sqrt(i) at the position i*rank in the vector vec_test, the other entries stay empty for each process
		vec_test[i] = sqrt(double(i));

#ifndef NO_MPI		
	MPI_Allreduce(MPI_IN_PLACE, &vec_test[0], size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  // this is the actual part, where the processes talk to each other. Here the vectors vec_test from all the processes are added and the sum is given back to every process. MPI_IN_PLACE ist the command to write the sum again in vec_test
#endif				
	
	if(myrank==0)														// use oinly the first process to print out stuff
	{
		for(int i=0; i<size; i++)
			cout << ", " << vec_test[i];
	}

#ifndef NO_MPI
	MPI_Finalize();
#endif	

	
}
