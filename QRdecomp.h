#include "nrutil.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cassert>
#include <vector>

using namespace std;
/*
routine of allocate matrix of m * n, 
if flag = true,row major layout
if flag = false, column major layout 
*/
vector<vector<float> > fmatrix(int m,int n,bool flag=true)
{
	int i;
	vector<vector<float> > matrix;
	int majorD,minorD;
	if(flag){
		majorD = m;
		minorD = n;
	}
	else{
		majorD = n;
		minorD = m;
	}
	for(i=0;i<majorD;i++){
		vector<float> row(minorD,0.0);
		matrix.push_back(row);
	}		

	return matrix;
}

/*
routine of allocate vector is 1...n
*/
vector<float> fvector(int n)
{
	vector<float> v(n,0.0);
	return v;
}

vector<vector<float> > eye(int m,int n)
{
	int i;
	vector<vector<float> > matrix = fmatrix(m,n,false);
	assert(m>n);
	for(i=0;i<n;i++)
		matrix[i][i] = 1;
	return matrix;
}

/*
Gram-Schmidt Orthogonalization process
*/
void orthogonal(vector<vector<float> > a,vector<float> c,int m,int n)
{
	int i,j,k;
	float sum = 0.0,scale,dotProduct;
	for(i=0;i<n;i++){
		for(j=0;j<m;j++)
			c[j] = a[i][j];
		for(j=i-1;j>=0;j--)
		{
			dotProduct = 0.0;
			for(k=0;k<m;k++)
				dotProduct+= c[k]*a[j][k];
			for(k=0;k<m;k++)
				a[i][k] -= (dotProduct * a[j][k]);
		}
		sum = 0.0;
		for(j=0;j<m;j++)
			sum += a[i][j]*a[i][j];
		scale = sqrtf(sum);
		for(j=0;j<m;j++)
			a[i][j] = a[i][j]/scale;
	}
	return;
}

/*
Constructs the QR decomposition of a[1..n][1..n]. 
The upper triangular matrix R is returned in the upper triangle of a, 
except for the diagonal elements of R which are returned in d[1..n]. 
The orthogonal matrix Q is represented as a product of n − 1 Householder matrices Q1...Qn−1,
where Qj =1−uj ⊗uj/cj. The ith component of uj is zero for i=1,...,j−1 
while the nonzero components are returned in a[i][j] for i = j, . . . , n. 
sing returns as true (1) if singularity is encountered during the decomposition, 
but the decomposition is still completed in this case; otherwise it returns false (0).
*/
void qrdcmp(float **a, int n, float *c,float *d,int* sing)
{
	int i,j,k;
	float scale,sigma,sum,tau;

	*sing = 0;
	for(k=1;k<n;k++){
		scale = 0.0;
		for(i=k;i<=n;i++) scale = FMAX(scale,fabs(a[i][k]));
		/*Singular case*/ 
		if(scale == 0.0) {
			*sing = 1;
			c[k] = d[k] = 0.0;

		} 
		/*Form Qk and Qk*A*/
		else
		{
			for(i=k;i<=n;i++) a[i][k]/= scale;
			for(sum = 0.0,i=k;i<=n;i++) sum+= SQR(a[i][k]);
			sigma = SIGN(sqrt(sum),a[k][k]);
			a[k][k] += sigma;
			c[k] = sigma*a[k][k];
			d[k] = -scale*sigma;
			for(j=k+1;j<=n;j++){
				for(sum=0.0,i=k;i<=n;i++) sum+= a[i][k]*a[i][j];
				tau = sum/c[k];
				for(i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
			}
		}
	}
	d[n] = a[n][n];
	if (d[n] == 0.0) *sing = 1;
}
