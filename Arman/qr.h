#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <vector>
#include <iostream>

/*
Gram-Schmidt Orthogonalization process for column-major matrices
*/

template <class T>
void orthogonal_col(std::vector<std::vector<T> > a, long nrl,long nrh,long ncl,long nch)
{
	std::vector<T> c,dummy;
    long i,j,k;
    
    T sum = 0.0,scale,dotProduct;
    
    /*for(j=0;j<=a.size();j++)
    {
    	std::copy(a[j].begin(), a[j].end(), std::ostream_iterator<T>(std::cout, " "));
  		std::cout << std::endl;
    }*/
    for(j=ncl;j<=nch;j++)
    {
        for(i=nrl;i<=nrh;i++)
            c.push_back(a[j][i]);
        std::cout << "vector " << j << "\n";
        std::copy(c.begin(), c.end(), std::ostream_iterator<T>(std::cout, " "));
  		std::cout << std::endl;  
        for(i=j-1;i>=0;i--)
        {
            dotProduct = 0.0;
            for(k=nrl;k<=nrh;k++)
                dotProduct+= c[k]*a[i][k];
            for(k=nrl;k<=nrh;k++)
                a[j][k] -= (dotProduct * a[i][k]);
        }
        sum = 0.0;
        for(i=nrl;i<=nrh;i++)
            sum += a[j][i]*a[j][i];
        std::cout << "Matrix on iteration " << j << " before scaling\n";
        std::copy(c.begin(), c.end(), std::ostream_iterator<T>(std::cout, " "));
  		std::cout << std::endl;  
        scale = sqrt(sum);
        std::cout << "Scale on iteration " << j << ": " << scale << "\n";
        for(i=nrl;i<=nrh;i++)
            a[j][i] = a[j][i]/scale;
        std::cout << "matrix on iteration " << j << "\n";
        std::copy(a[j].begin(), a[j].end(), std::ostream_iterator<T>(std::cout, " "));
  		std::cout << std::endl;  
        c.clear();
        c.shrink_to_fit();
    }
    
    
    return;
}


/*
Gram-Schmidt Orthogonalization process
*/

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
