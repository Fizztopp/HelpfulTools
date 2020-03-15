// Uses FFTW to calculate the real+imaginary part of a test function Testf() by input of its real/imaginary part <-- the amount of information doesn't change!!
// (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

//g++ -O2 -std=c++11 FFT_sample.cc -lfftw3 -lm


#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <fftw3.h>

#define N_FFTW 2000														// size of input/output arrays for FFTW
#define Omega_min -200.                                                  // min frequency for FFTW
#define Omega_max +200.                                                  // maximum frequency for FFTW


using namespace std;

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);


cdouble Testf(double omega, double eps, double eta)
/*
 test function
 */
{
	return 1./(omega - eps + II*eta);
}	

//IMPORTANT: No normalization in FFTW -> foward + backward transform multiplies factor N_FFTW

int main()
{
	double eps = 5.0;
	double eta = 0.5;

	fftw_complex input[N_FFTW];											// input array
	
	fftw_plan pf;                                                       // forward plan
	fftw_plan pb;														// backward plan
	
	pf = fftw_plan_dft_1d(N_FFTW, input, input, FFTW_FORWARD, FFTW_MEASURE);
	pb = fftw_plan_dft_1d(N_FFTW, input, input, FFTW_BACKWARD, FFTW_MEASURE);
	
	
	dvec OMEGA(N_FFTW);													// array to store frquencies
	
	
	for(int iw=0; iw<N_FFTW; iw++)
	{
		OMEGA[iw] = Omega_min + iw*(Omega_max-Omega_min)/N_FFTW;        // initialize Omega		
		input[iw][0] = Testf(OMEGA[iw], eps, eta).real();		        // put f(Omega).real() as input array ([0] indicates real part for fftw_complex)
	}
		
	fftw_execute(pf);													// FT[Testf(t)]
	
	for(int it=N_FFTW/2; it<N_FFTW; it++)
	{	
		input[it][0] = 0;       										// deleta one half of the array	<- Theta(t)
		input[it][1] = 0;      															
	}

	
	fftw_execute(pb);													// FT^(-1)[ FT[Testf(t)] * Theta(t)]
		
	fftw_destroy_plan(pf);	
	fftw_destroy_plan(pb);
	
	cvec RESULT(N_FFTW);
	
	ofstream myfile ("test_data.txt");
	 
	for(int iw=0; iw<N_FFTW; iw++)
	{
		input[iw][0] *= 2./(double(N_FFTW));                            // factor 2: The FT of a real function (imaginary part of Testf()) will be symmetric, but since we delete the half pf the array we have to multiply by 2
		input[iw][1] *= 2./(double(N_FFTW));                            
		//cout << input[it][0] << " " << input[it][1] << endl;
		//cout << Testf(OMEGA[iw], eps, eta).real() << " " << Testf(OMEGA[iw], eps, eta).imag() << endl << endl;
		myfile << input[iw][0] << " " <<  input[iw][1];
		myfile << endl;
	}
	myfile.close();
}	
