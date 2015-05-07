#ifndef SPARSE_MATRIX_H
#define	SPARSE_MATRIX_H

#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>

template <class T>
class Matrix
{
public:
    typedef std::map<size_t, std::map<size_t , T> > mat_t;
    typedef typename mat_t::iterator row_iter;
    typedef std::map<size_t, T> col_t;
    typedef typename col_t::iterator col_iter;

    Matrix(size_t i) { m=i; n=i; } //constructor of diagonal matrix
    Matrix(size_t i, size_t j) { m=i; n=j; }  //constructor of general matrix
	
    inline
    T& operator()(size_t i, size_t j)
    {
        if(i>=m || j>=n) throw;
        return mat[i][j];
    }
    inline
    T operator()(size_t i, size_t j) const
    {
        if(i>=m || j>=n) throw;
        return mat[i][j];
    }
    
    const mat_t Map()
    {
    	mat_t mat_out;
    	mat_out = this->mat;
    	return mat_out;
    }
    
    int Cols(){return n;}
    int Rows(){return m;}
    
    void Set_Cols(int i){n=i;}
    void Set_Rows(int i){m=i;}
    
    void operator=(const Matrix& other)
	{
		mat = other.mat;
		n = other.n;
		m = other.m;
	}

    std::vector<T> operator*(const std::vector<T>& x)
    {  //Computes y=A*x
        if(this->m != x.size()) throw;

        std::vector<T> y(this->m);
        T sum;

        row_iter ii;
        col_iter jj;

        for(ii=this->mat.begin(); ii!=this->mat.end(); ii++)
        {
            sum=0;
            for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
                sum += (*jj).second * x[(*jj).first];
            y[(*ii).first]=sum;
        }

        return y;
    }
	
	std::vector<T> Coopmans_guess()
    {  //Extracts diagonal of the matrix

        std::vector<T> y, z;
		
		T dummy;
		int i, j;

        row_iter ii;
        col_iter jj;
        
        //std::cout << "Forming Coopmans guess, size of Hamiltonian: "
        // << this->m << std::endl;
		
        for(ii=this->mat.begin(); ii!=this->mat.end(); ii++)
            for(jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
                if((*ii).first==(*jj).first)     //debug it!!
                {
					y.push_back((*jj).second);
					//std::cout << (*ii).first << "\t" 
                	//<<(*jj).first <<"\t" << y.back() << std::endl;
				}
		//std::cout << "Size of diagonal of Hamiltonian is: " << y.size() << std::endl;
        	
		return y;
    }
    
    void printMat()
    {
        row_iter ii;
        col_iter jj;
        for(ii=this->mat.begin(); ii!=this->mat.end(); ii++){
            for( jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++){
                std::cout << (*ii).first << ' ';
                std::cout << (*jj).first << ' ';
                std::cout << (*jj).second << ' ';
            }
            std::cout << std::endl;
        } std::cout << std::endl;
    }
    
    /*std::ostream &MatrixPrint(std::ostream &out) 
	{
		row_iter ii;
    	col_iter jj;
    	for(ii=this->mat.begin(); ii!=this->mat.end(); ii++)
    	{
         	for( jj=(*ii).second.begin(); jj!=(*ii).second.end(); jj++)
         	{
            	out << (*ii).first << ' ';
              	out << (*jj).first << ' ';
              	out << (*jj).second << std::endl;
        	}
    	} 
    	out << std::endl;
    	return out;
	}
    
	friend std::ostream &operator<<(std::ostream &out, const Matrix &mtrx) 
	{
		out = mtrx.MatrixPrint(out);
    	return out;
	}*/
    
    Matrix(){}

private:
    mat_t mat;
    size_t m;     //row
    size_t n;     //column
};

#endif	/* _SPARSE_MATRIX_H */