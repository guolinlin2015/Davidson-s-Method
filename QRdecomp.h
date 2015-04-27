#ifndef QR_H
#define QR_H

#include "nrutil.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cassert>


float **eye(long nrl,long nrh,long ncl, long nch)
{
	long i;
	float **m = matrix(nrl,nrh,ncl,nch);
	long nrow = nrh-nrl +1,ncol = nch-ncl +1;
	assert(nrow>=ncol);
	for(i=ncl;i<=nch;i++)
		m[i][i] = 1.0;
	return m;
}

/*
Gram-Schmidt Orthogonalization process
*/
void orthogonal(float **a,float *c,long nrl,long nrh,long ncl,long nch)
{
	long i,j,k;
	float sum = 0.0,scale,dotProduct;
	for(j=ncl;j<=nch;j++){
		for(i=nrl;i<=nrh;i++)
			c[i] = a[i][j];
		for(i=j-1;i>=ncl;i--)
		{
			dotProduct = 0.0;
			for(k=nrl;k<=nrh;k++)
				dotProduct+= c[k]*a[k][i];
			for(k=nrl;k<=nrh;k++)
				a[k][j] -= (dotProduct * a[k][i]);
		}
		sum = 0.0;
		for(i=nrl;i<=nrh;i++)
			sum += a[i][j]*a[i][j];
		scale = sqrtf(sum);
		for(i=nrl;i<=nrh;i++)
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

#endif
