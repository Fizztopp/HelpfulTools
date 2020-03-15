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

#ifndef NO_OMP
    #include <omp.h>
#endif

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);

int main (){
	//************** OPEN_MP INIT **************************************
#ifndef NO_OMP 	  
	cout << "# of processes " << omp_get_num_procs() << endl;
#pragma omp parallel 
	cout << "Thread " << omp_get_thread_num() << " out of " << omp_get_num_threads() << " says hello!" << endl;     
#endif
	//******************************************************************
	
	double dtime, dtime0, dtime1;
	
	// create random Hermirian matrix	
	dvec *TEMP0 = new dvec(2);
	dvec *TEMP1 = new dvec(2);
	dvec *TEMP2 = new dvec(2);
	dvec *TEMP3 = new dvec(2);

	dvec *p1, *p2, *p3;

	(*TEMP1)[0] = 1.0;
	(*TEMP1)[1] = 2.0;
	
	(*TEMP2)[0] = 3.0;
	(*TEMP2)[1] = 4.0;

	(*TEMP3)[0] = 5.0;
	(*TEMP3)[1] = 6.0;

	p1 = TEMP1;
	p2 = TEMP2;
	p3 = TEMP3;
	
	TEMP1 = p2;
	TEMP2 = p3;
	TEMP3 = p1;
	
	cout << endl << "Pointers get same address" << endl;
	cout << "TEMP1[0]: " << (*TEMP1)[0] << " " << "TEMP1[1]: " << (*TEMP1)[1] << endl;
	cout << endl << "TEMP2[0]: " << (*TEMP2)[0] << " " << "TEMP2[1]: " << (*TEMP2)[1] << endl;
	cout << endl << "TEMP3[0]: " << (*TEMP3)[0] << " " << "TEMP3[1]: " << (*TEMP3)[1] << endl;
	
	(*TEMP3)[0] = 50.0;
	(*TEMP3)[1] = 60.0;
	
	cout << endl << "Pointers get same address" << endl;
	cout << "TEMP1[0]: " << (*TEMP1)[0] << " " << "TEMP1[1]: " << (*TEMP1)[1] << endl;
	cout << endl << "TEMP2[0]: " << (*TEMP2)[0] << " " << "TEMP2[1]: " << (*TEMP2)[1] << endl;
	cout << endl << "TEMP3[0]: " << (*TEMP3)[0] << " " << "TEMP3[1]: " << (*TEMP3)[1] << endl;
	
	dvec *TEMP4 = new dvec(2);
	dvec *TEMP5 = new dvec(2);
	
	(*TEMP4)[0] = 1.0;
	(*TEMP4)[1] = 2.0;

	TEMP5[0] = TEMP4[0];
	
	(*TEMP4)[0] = 3.0;
	(*TEMP4)[1] = 4.0;
	
	cout << endl << "Pointers get same values" << endl;
	cout << "TEMP4[0]: " << (*TEMP4)[0] << " " << "TEMP4[1]: " << (*TEMP4)[1] << endl;
	cout << endl << "TEMP5[0]: " << (*TEMP5)[0] << " " << "TEMP5[1]: " << (*TEMP5)[1] << endl;
	double a = 0.1;
	double *ap = &a; 
	double *bp;
	bp = ap; 
	a = 0.2;

	cout << endl << "a: " << a << " " << "ap: " << *ap << " " << "bp: " << *bp << endl;
	
	delete TEMP0, TEMP1, TEMP2, TEMP3, TEMP4, TEMP5;
}

