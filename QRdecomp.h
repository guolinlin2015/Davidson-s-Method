#ifndef QR_H
#define QR_H

#include "mynrutil.h"
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

#endif
