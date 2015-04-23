#include <iostream>
#include <ctime>
#include <cstdlib>
#include "QRdecomp.h"
#include "mvutil.h"

using namespace std;

int main(int argc,const char* argv[])
{
	int n = 500;
	float tol = 1e-8;
	int k = 4;
	int l = 2*k;
	int i,j;
	float sparsity = 0.0000001;

	srand((unsigned)time(NULL));

	//Construct Matrix
	// A & AT: row major
	vector<vector<float> > A = fmatrix(n, n);
	for(i=0;i<n;i++){
		A[i][i] = i+1;
	}
	for(i=0;i<n;i++){
		for(j=0;j<n;j++)
			A[i][j] = A[i][j] + sparsity * rand();
	}
	vector<vector<float> > AT = fmatrix(n, n);
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			AT[i][j] = A[j][i];

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			A[i][j] = (A[i][j]+AT[i][j])/2;

	//Construct Guess Vector
	//B : column major
	vector<vector<float> > B = eye(n, l);
	vector<float> c = fvector(n);
	orthogonal(B,c,n,l);

	//Subspace matrix
	//Abar : row major
	vector<vector<float> > Abar = fmatrix(l,l);
	for(i=0;i<l;i++){
		for(j=0;j<l;j++){
			Abar[i][j] = dotProduct(B[i],MVM(A,B[j],n,l),n);
		}
	}







}