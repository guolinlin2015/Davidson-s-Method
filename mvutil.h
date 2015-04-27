#ifndef MVUTIL_H
#define MVUTIL_H

#include "mynrutil.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

float dotProduct(float *a,float *b,long nl,long nh){
	float sum = 0.0;
	for(long i=nl;i<=nh;i++){
		sum += (a[i]*b[i]);
	}
	return sum;
}

float *VectorAdd(float *a,float *b,long nl,long nh){
	float *result = vector(nl, nh);
	for (long i=nl;i<=nh;i++)
		result[i] = a[i] + b[i];
	return result;
}

float *VectorSubtract(float *a,float *b,long nl,long nh){
	float *result = vector(nl, nh);
	for(long i=nl;i<=nh;i++){
		result[i] = a[i] - b[i];
	}
	return result;
}

float *ScaleVector(float *a,float scale,long nl,long nh){
	float *result = vector(nl, nh);
	for(long i = nl;i<=nh;i++){
		result[i] = a[i]*scale;
	}
	return result;
}

float VectorNorm(float *a,long nl,long nh){
	float sum = 0.0;
	for (long i=nl;i<=nh;i++){
		sum += a[i];
	}
	return sqrtf(sum);
}

void VectorNormalize(float *a,long nl,long nh){
	float norm = VectorNorm(a, nl, nh);
	for(long i=nl;i<=nh;i++)
		a[i] = a[i]/norm;
	return;
}

float *MatrixGetRow(float **m,long nrl,long nrh,long nc){
	float *v = vector(nrl, nrh);
	for(long i = nrl;i<=nrh;i++)
		v[i] = m[i][nc];
	return v;
}

float **MatrixSubtract(float **a,float **b,long nrl,long nrh,long ncl,long nch){
	float **result = matrix(nrl, nrh, ncl, nch);
	for(long i=nrl;i<=nrh;i++)
		for(long j=ncl;j<=nch;j++)
			result[i][j] = a[i][j] - b[i][j];

	return result;
}

float **VectorVectorTM(float *a,float *b,long nrl,long nrh,long ncl,long nch){
	float **result = matrix(nrl, nrh, ncl, nch);
	for(long i=nrl;i<=nrh;i++)
		for(long j=ncl;j<=nch;j++)
			result[i][j] = a[i] * b[j];
	return result;
}

float **MatrixMatrixM(float **a,float **b,long nrl,long nrh,long ncl,long nch){
	float **result =  matrix(nrl, nrh, ncl, nch);
	for(long i=nrl;i<=nrh;i++)
		for(long j=ncl;j<=nch;j++){
			float sum = 0.0;
			for(long k=ncl;k<=nch;k++){
				sum += a[i][k] * b[k][j];
			}
			result[i][j] = sum;
		}
	return result;
}

float **MatrixAddCol(float **m,float *v,long nrl,long nrh,long ncl,long nch){
	float **result = matrix(nrl,nrh,ncl,nch+1);
	for(long j=ncl;j<=nch;j++)
		for(long i=nrl;i<=nrh;i++)
			result[i][j] = m[i][j];

	for(long i=nrl;i<=nrh;i++)
		result[i][nch+1] = v[i];
	return result;
}

//a must be row major
float *MVM(float **a,float *b,long nrl,long nrh,long ncl,long nch){
	float *v = vector(nrl,nrh);
	memset(v,0,sizeof(float)*(nrh-nrl+1));
	for(long i=nrl;i<=nrh;i++){
		v[i] = dotProduct(a[i],b,ncl,nch);
	}
	return v;
}

#endif