#ifndef MVUTIL_H
#define MVUTIL_H

#include <vector>
#include <iostream>

using namespace std;


float dotProduct(vector<float> a,vector<float> b,int n){
	float sum = 0.0;
	for(int i=0;i<n;i++){
		sum += (a[i]*b[i]);
	}
	return sum;
}

//a must be row major
vector<float> MVM(vector<vector<float> > a,vector<float> b,int m,int n){
	vector<float> v(m,0.0);
	for(int i=0;i<m;i++){
		v[i] = dotProduct(a[i],b,n);
	}
	return v;
}

#endif