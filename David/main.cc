#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>

#include "davidson.h"
#include "vector_algebra.h"
#include "sparse_matrix.h"
#include "eigen.h"

#define NInitGuess 8
#define N_eig 4

using namespace std;

/*template <class T>
Matrix<T> init(Matrix<T> &A)
{
    ifstream fin("bcsstk11.mtx");
    ofstream fout("bcsstk11.out");
    cout << "\t Opening file...";
    int m,n,length,i,j;
    if(fin.is_open())
    {
    	cout << "\tSuccess\n";
    	
        fin >> n;
        fin >> m;
		fin >> length;
		
		A.Set_Rows(n);
		A.Set_Cols(m);
		
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
    
    return A;
}*/

template <class T>
Matrix<T> init(Matrix<T> &A)
{
	int i,j;
	
	int n=50;
	int length=n*n;
	double sparsity = 0.95;
	
	A.Set_Rows(n);
	A.Set_Cols(n);
	
	srand (time(NULL));
	
	for(i=0; i<n; i++)
	{
		for(j=0; j<i; j++)
        	if( rand()/(double)RAND_MAX > ( 1-(1-sparsity)/2 ) )
        	{
        		A(i,j)=rand()/(double)RAND_MAX;
        		A(j,i)=A(i,j);
        	}
        A(i,i)=5*i+10;
    }
    
    return A;
}

int main()
{

    Matrix<long double> A;
    
    long double sum;
    ofstream output;
    int n,k=0;
    clock_t t0,tf;
    cout.precision(17);
    
    cout << "Initializing matrix...\n";
    
    A=init(A);
    
    vector<long double> x(A.Cols(),1), y;
    
    cout << "Success\n";
    
    A.printMat();
    
    Davidson_parameters params(4, 8, 40, 4, "Coopmans") ;
    
    output.open("davidson.out");
    
	Davidson<long double> davidson(A, params);
	
	davidson.solve(output);
	
		
	/*//timed multiplication
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
    cout << "\n" << tmp <<"\n";*/

    return 0;
}