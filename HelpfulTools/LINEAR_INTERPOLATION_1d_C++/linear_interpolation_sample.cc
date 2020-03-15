// Routine for 1d linear interpolation
// (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de
//g++ -O2 -std=c++11 tests.cc -lfftw3 -lm

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include <math.h>
#include <fftw3.h>

#define N 1000														// size of input/output arrays for FFTW
#define x_min -5.                                                // min frequency for FFTW
#define x_max +5.                                                // maximum frequency for FFTW
#define x_test 2.12345

using namespace std;
typedef vector<double> dvec;                     					    // vectors with real double values

double Testf(double x)
{
	return x*x/2.;
}	

inline double linear_interpolation(dvec &fx, dvec &x, double x_target) 
{     
/*
 -fx: vector of function values
 -x: vector of function arguments
 -x_target: target argument
 */	   
    vector<double>::iterator iter0; 
    vector<double>::iterator iter1;
    int i0;
    int i1;
	iter0 = lower_bound(x.begin(), x.end(), x_target);
	i0 = (iter0 - x.begin())-1;
	i1 = i0 + 1;
    return fx[i0] + (fx[i1] - fx[i0])/(x[i1] - x[i0])*(x_target - x[i0]); 
} 

int main()
{
	dvec X(N);
	dvec FX(N);													
	
	for(int iw=0; iw<N; iw++)
	{
		X[iw] = x_min + iw*(x_max-x_min)/N;        
		FX[iw] = Testf(X[iw]);
		//cout << X[iw] << " " << FX[iw] << endl;
	}
	cout << "Interpolation value: " << linear_interpolation(FX, X, x_test) << endl;
	cout << "Analytic value: " << Testf(x_test) << endl; 
	cout << "Absolute deviation: " << abs(linear_interpolation(FX, X, x_test)-Testf(x_test) ) << endl;
}	
