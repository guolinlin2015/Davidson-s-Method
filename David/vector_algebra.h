#ifndef VECTOR_ALGEBRA_H
#define	VECTOR_ALGEBRA_H

#include <cstdlib>
#include <map>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator

template <class T>
std::vector<std::vector<T> > matrix_transpose (const std::vector<std::vector<T> > &A) 
{
  	std::vector<std::vector<T> > result;
  	std::vector<T> row;
  	int i,j;
  
  	for(i=0; i<A[0].size(); i++)
  	{
  		for(j=0; j<A.size(); j++)
  			row.push_back(A[i][j]);
  		result.push_back(row);
  		//std::copy(row.begin(), row.end(), std::ostream_iterator<T>(std::cout, " "));
  		//std::cout << std::endl;
  		row.clear();
  		row.shrink_to_fit();
  		
  	}
  	
  	
  	
  	return (result);
}

template <class T>
T dot (const std::vector<T> &a, const std::vector<T> &b) 
{
  	T result;
  	result = 0;
  	int i;
  	
  	if(a.size()!=b.size())
  	{
  		std::cout << "Vector-vector dot product mismatch!\n";
  		throw;
  	}
  
  	for(i=0; i<a.size(); i++)
  		result+=a[i]*b[i];
  	return (result);
}

template <class T>
std::vector<T> vector_inverse (
const std::vector<T> &A) 
{
  	std::vector<T> result;
  	int i;
  
  	for(i=0; i<A.size(); i++)
  	{	
  		if(A[i]!=0)
  			result.push_back(1/A[i]);
  		else
  		{
  			std::cout << "got 0 in the preconditioner inverse!\n";
  			throw;
  		}
  	}
  	return (result);
}
template <class T>
std::vector<T> diag_vec_mult(
std::vector<T> a,std::vector<T> b)  
//a is vector of digonal elements of some matrix, b is some vector
{
	std::vector<T> result;
  	int i;
    
    /*if(a.size()!=b.size())
  	{
  		std::cout << "Vector-vector 'strange multiplication' dimentions mismatch!\n";
  		std::cout << "Vector width: " << a.size() <<"\n";
  		std::cout << "Vector size: " << b.size() <<"\n";
  		throw;
  	}*/
    
  	for(i=0; i<a.size(); i++)
  		result.push_back(a[i]*b[i]);
  		
  	return (result);
}

template <class T>
T norma (const std::vector<T> &a) 
{
  	T result;
  	result=sqrt(dot(a,a));
  		
  	return (result);
}

template <class T>
std::vector<std::vector<T> > operator*( 
const std::vector<std::vector<T> > &A, const T &b) 
{
	std::vector<std::vector<T> > result;
  	int i, j;
  
  	for(i=0; i<A.size(); i++)
  	{
  		result.push_back(A[i]);
  		for(j=0; j<A[0].size(); j++)
  			result[i][j]=A[i][j]*b;
  	}	
  	return (result);
}

template <class T>
std::vector<T> operator*( 
const std::vector<std::vector<T> > &A, const std::vector<T> &v) 
{
	std::vector<T> result;
  	int i,j;
  	
  	if(A.size()!=v.size())
  	{
  		std::cout << "Matrix-vector multiplication dimentions mismatch!\n";
  		std::cout << "Matrix width: " << A[0].size() <<"\n";
  		std::cout << "Matrix height: " << A.size() <<"\n";
  		std::cout << "Vector size: " << v.size() <<"\n";
  		throw;
  	}
  
  	for(i=0; i<A[0].size(); i++)
  	{
  		result.push_back(0);
  		for(j=0; j<v.size(); j++)
  			result[i]+=A[j][i]*v[j];
  	}		
  	//std::cout << "Resulting vector size: " << result.size() <<"\n";
  	
  	return (result);
}

template <class T>
T operator*( 
const std::vector<T> &a, const std::vector<T> &b) 
{
	T result = 0;
  	int i;
  
  	for(i=0; i<a.size(); i++)
  		result+=(a[i]*b[i]);
  	return result;
}

template <class T>
std::vector<T> operator*( 
const std::vector<T> &v, const T A) 
{
	std::vector<T> result;
  	int i,j;
    
    
  	for(i=0; i<v.size(); i++)
  		result.push_back(A*v[i]);
  		
  	return (result);
}

template <class T>
std::vector<T> operator/( 
const std::vector<T> &v, const T A) 
{
	std::vector<T> result;
  	int i,j;
  
  	for(i=0; i<v.size(); i++)
  		result.push_back(v[i]/A);
  	return (result);
}

template <class T>
std::vector<std::vector<T> > operator-( 
const std::vector<std::vector<T> > &A, const T &b) 
{
	std::vector<std::vector<T> > result;
  	int i;
  
  	for(i=0; i<A.size(); i++)
  	{
  		result.push_back(A[i]);
  		result[i][i]=A[i][i]-b;
  	}	
  	return (result);
}

template <class T>
std::vector<T> operator-( 
const std::vector<T> &a, const std::vector<T> &b) 
{
	std::vector<T> result;
  	int i;
  
  	for(i=0; i<a.size(); i++)
  		result.push_back(a[i]-b[i]);
  		
  	return (result);
}

template <class T>
std::vector<T> operator-( 
const std::vector<T> &a, const T &b) 
{
	std::vector<T> result;
  	int i;
  
  	for(i=0; i<a.size(); i++)
  		result.push_back(a[i]-b);
  		
  	return (result);
  	
}


/*
template <class T>
bool operator>(const std::vector<T> &v, T &tol) 
{
	return ( sqrt( dot(v,v) ) > tol );
}

template <class T>
bool operator<(const std::vector<T> &v, T &tol) 
{
	return ( sqrt( dot(v,v) ) < tol);
}*/
#endif	/* VECTOR_ALGEBRA_H */