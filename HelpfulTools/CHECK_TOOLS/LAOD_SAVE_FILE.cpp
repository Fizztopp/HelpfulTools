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

typedef complex<double> cdouble;                  						// typedef existing_type new_type_name ;
typedef vector<double> dvec;                     					    // vectors with real double values
typedef vector<cdouble> cvec;                     						// vectors with complex double values

cdouble II(0,1);

int main (){
	cvec TEST(6);
	
	ifstream in("TEST_ARRAY.dat");
	if (!in) 
	{
		cout << "Cannot open file.\n";
		return 0;
	}
	for(int i=0; i<6; i++)
		in >> TEST[i];			
		
	in.close();
	in.clear();
	
	ofstream myfile("TEST_ARRAY_1.dat");
	if (myfile.is_open())
	{
		for(int i=0; i<6; i++){
			myfile << TEST[i] << " ";
		}	
		myfile.close();	
	}
	else cout << "Unable to open file" << endl;	
}

