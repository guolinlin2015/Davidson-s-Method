#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include "sparse_matrix.h"

using namespace std;

int init(Matrix<long double> &A)
{
    ifstream fin("bcsstk11.mtx");
    cout << "\t Opening file...";
    int n,m,length,i,j;
    if(fin.is_open())
    {
    	cout << "\tSuccess\n";
    	
        fin >> n;
        fin >> m;
		fin >> length;
		
		cout << n << "\t" << m << "\t" << length << "\n";
		
        int k=0;
        while(k<length)
        { 
        	fin >> i; 
        	fin >> j; 
        	fin >> A(i-1,j-1); 
        	//cout << i << "\t" << j  << "\t" << A(i-1,j-1) << "\n";
        	k++; 
        }
        fin.close();
    }
    else return 0;
    return n;
}

int main()
{
    Matrix<long double> A(1806);
    vector<long double> x(1806,1), y;
    long double sum;
    int n,k=0;
    clock_t t0,tf;
    cout.precision(17);
    
    cout << "Initializing matrix...\n";
    
    n=init(A);
    
    cout << "Success\n";
    
    //timed multiplication
    t0=clock();
    for(int i=0;i<1000;i++){ y=A*x; }
    tf=clock()-t0;
    
    //output last element
    cout << y[1805];
    //output time to multiply
    cout << "\n\n" << (double)tf/((double)(CLOCKS_PER_SEC)*(1000.0));

    //do naive multiply for last row to check last element of y
    long double tmp=0;
    for(int i=0;i<1806;i++) tmp+=A(1805,i);
    cout << "\n" << tmp <<"\n";

    return 0;
}