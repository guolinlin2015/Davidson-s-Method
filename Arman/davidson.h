#ifndef DAVIDSON_H
#define	DAVIDSON_H

#include <cstdlib>
#include <map>
#include <vector>
#include <iostream>
#include "vector_algebra.h"
#include "sparse_matrix.h"
#include "eigen.h"
#include "qr.h"

using namespace std;

template <class T>
class Davidson;

class Davidson_parameters
{
public:
	Davidson_parameters(
	size_t eig, size_t nguess, size_t nvecmax, 
	int tolerance, string guess) 
	{
		N_eig=eig; N_guess=nguess; Tol=tolerance; 
		Nvec_max=nvecmax; Guess_type=guess; 
	}
	
	void set(const Davidson_parameters& other)
	{
		N_eig=other.N_eig; N_guess=other.N_guess; 
		Tol=other.Tol; Nvec_max=other.Nvec_max; 
		Guess_type=other.Guess_type;
	}
	
	void print()
	{
		cout << N_eig << " is Number of requested eigen states\n";
		cout << N_guess << " is Number of initial guess vectors \n";
		cout << Nvec_max << " is Maximum number of projected space \n";
		cout << Tol << " is Tolerance for Davidson convergence \n";
		cout << Guess_type << " is Guess type \n";
	}
	
	Davidson_parameters(){}
	
	template <class T>
	friend class Davidson;

private:
	size_t N_eig;
	size_t N_guess;
	string Guess_type;
	size_t Nvec_max;
	int Tol;
};


template <class T>
class Davidson
{
public:
	
	//Davidson inititializing

    Davidson(const Matrix<T> &Hamiltonian, 
    const Davidson_parameters &parameters)
    {
	int i = 0, j = 0, k = 0;
	
	T dummy;
	
	cout << "Construction started \n";
	
    Hmltn = Hamiltonian;
    params.set(parameters);
    
	cout << "Matrix and parameters passed \n";
	
	params.print();
    
	sub_space_size=params.N_guess;
	converged = false;
	
		
	for(i=0; i<sub_space_size; i++)
		lambda.push_back(0.0);
	
	cout << "Lambda vector initialized \n";
	
	cout << "Diagonal initialized \n";
	
	diag = Hmltn.Coopmans_guess();	
	
	cout << "Diagonal is: \n";
	
	print_vec(diag,cout);
	
	vector<T> y;
	
	y=diag;
	
	vector<int> lowest(Hmltn.Cols());
	
	for(i=0; i<lowest.size(); i++)
		lowest[i]=i;
	
	for(i=0; i<y.size(); i++)
			for (j=i+1; j<y.size(); j++)  //sort Coopmans guess
				if (y[i]>y[j])
				{
					//std::cout << y[i] << " " << y[j] << " \n";
					dummy = y[i];
					y[i] = y[j];
					y[j] = dummy;
					k = lowest[i];
					lowest[i] = lowest[j];
					lowest[j] = k;	
				}
	
	v.push_back(lambda);
	for(i=sub_space_size; i<Hmltn.Cols(); i++)
		v[0].push_back(0.0);
	
	for(i=1; i<sub_space_size; i++)
		v.push_back( v[0] );
	
	for(i=0; i<sub_space_size; i++)
		v[lowest[i]][lowest[i]]=1.0;	
	
	
	cout << "Guess vectors initialized \n";	
	
	diag = Hmltn.Coopmans_guess();        //form diagonal of Hamiltonian
	
	sigma_v=v;
	
	cout << "Sigma vectors initialized \n";
	
	for(i=0; i<params.N_eig; i++)
	{
		residual.push_back(v[0]);
		eig_vec.push_back(lambda);  //eigen vec and residuals
	}
	
	cout << "Residuals and eigen vectors initialized \n";
	
	for(i=0; i<sub_space_size; i++)
		A_proj.push_back( lambda );
	
	cout << "Projected Hamiltonian initialized \n";
	
	
}
    
    void solve(ofstream& output)
    {
	
		initialize();
	
		while(converged!=true)
			next_step(output);
	
		finalize(output);
	}
	
    Davidson(){}
    
    ~Davidson(){}
    
private:

    /**    \brief Initializes the solver
     **/
    void initialize()
	{	
	int i = 0, j = 0;
	
	cout << "Initialization started \n";
	
	for(i=0; i<sub_space_size; i++)
		sigma_v[i]=Hmltn*v[i];   // sigma_v
		
	cout << "Sigma vectors formed \n";
		
	for(i=0; i<sub_space_size; i++)
		for(j=0; j<sub_space_size; j++)
			A_proj[i][j] = dot(v[i],sigma_v[j]);   // A_proj
	
	cout << "Projected Hamiltonian formed \n";
			
	
	solve_subspace();   //updating eig_vec and lambda
	
	print_current_step(cout);
	
	}

    /**    \brief Performs one iteration step
     **/
    void next_step(ofstream& output)
	{
		int i, j;
	
		cout << "Step of Davidson \n";
	
		compute_residuals();
	
		print_current_step(output);
	
		if(converged==true)
			return;
	
		if(sub_space_size>params.Nvec_max)
		{
			collapse_subspace();
			return;
		}
	
		update_vectors();
	
	}


    /**    \brief Finishes the work
     **/
     
    void finalize(ofstream& output)
    {
	output << "Davidson converged \n";
	print_vec(lambda, output);
	output.close();
	return;
	}


    /**    \brief Collapses the previous step's subspace 
    		to the minimum number of vectors
     **/
     
    void collapse_subspace()
    {
		int i;
		sub_space_size = params.N_eig;
		v.clear();

		for(i=0; i<params.N_eig; i++)
			v.push_back(eig_vec[i]);
	}
    
	void compute_residuals()
	{
		int i=0, cnvg=0;
	
		cout << "Computing residuals\n";
	
		cout << "Residuals height: " << residual.size() <<"\n";
	 	cout << "Residuals width: " << residual[i].size() <<"\n";
	 	cout << "sigma_v height: " << sigma_v.size() <<"\n";
	 	cout << "sigma_v width: " << sigma_v[i].size() <<"\n";
	 	cout << "v height: " << v.size() <<"\n";
	 	cout << "v width: " << v[i].size() <<"\n";
	 	cout << "eig_vec height: " << eig_vec.size() <<"\n";
	 	cout << "eig_vec width: " << eig_vec[i].size() <<"\n";
	 	
		for(i=0; i<params.N_eig; i++)
		{
			residual[i]=sigma_v*eig_vec[i]-(v*eig_vec[i])*lambda[i];
			print_vec(residual[i], cout);
			if (norma(residual[i])>10.0^(-params.Tol))
				cnvg++;
		}
		
		if (cnvg==0)
			converged=true;
		
		cout << "Residuals calculated, converged:" << params.N_eig-cnvg << "\n";
	}

	
    void update_vectors()
    {
	std::vector<std::vector<T> > precon_res;
	std::vector<std::vector<T> > trial;
	int i, j=0, cnvg=0;
	
	cout << "Updating vectors\n";
	
	
	
	for(i=0; i<params.N_eig; i++)
		if (norma(residual[i])>pow(10.0,-params.Tol))
		{
			precon_res.push_back(
			diag_vec_mult(vector_inverse(diag-lambda[i]),residual[i]));  
			//diag is vector, lambda i is number, residual i is vector, 
			//product of multiplication must be vector!!
			cout << "Norm of residual:" << norma(residual[i])  << "\n";
		}
	cout << "Preconditioner applied to the residuals\n";
	
	trial = v;
	for(i=0; i<precon_res.size(); i++)
		trial.push_back(precon_res[i]);
	
	print_matrix(trial,cout);
	
	cout << "trial[0].size(): " << trial[0].size()  << "\n";
	cout << "trial.size(): " << trial.size()  << "\n";
	
	orthogonal_col(trial, 0, trial[0].size(), sub_space_size, trial.size() );
	
	print_matrix(trial,cout);
	
	cout << "Orthogonalized preconditioned residuals\n";
	
	for(i=0; i<precon_res.size(); i++)
	{
		if (norma(trial[i+sub_space_size])>10.0^(-params.Tol))
		{
			v.push_back(trial[i+sub_space_size]);
			print_vec(v.back(), cout);
			sigma_v.push_back(Hmltn*v.back());
			print_vec(sigma_v.back(), cout);
		}
	}
	
	cout << "Updated v and sigma_v\n";
	
	for(i=0; i<sub_space_size; i++)
		for(j=0; j<precon_res.size(); j++)
			A_proj[i].push_back(v[i]*(sigma_v[j]));   
			//update old rows (adding new cols)
	
	for(i=sub_space_size; i<sub_space_size+precon_res.size(); i++)
	{
		vector<T> row;
		for(j=0; j<sub_space_size+precon_res.size(); j++)
			row.push_back(v[j]*sigma_v[i]);   //add new elongated rows
		print_vec(row, cout);
		A_proj.push_back(row);   
	}
	cout << "Updated A_proj\n";
	
	sub_space_size+=precon_res.size();
	
	cout << "A_proj: dimentions: " << A_proj.size() << " by (";
	for(i=0; i< A_proj.size(); i++) 
		cout << A_proj[i].size() << " ";
	cout << ")\n";
	
	print_matrix(A_proj, cout);
	
	throw;
	
	solve_subspace();

	precon_res.clear();
	}
	
    void solve_subspace()
    {
	int i, j, k;
	double dummy;
	double **h; /* Hamiltonian matrix */
	double *d; /* Eigenvalues */ 
	double *e; /* Work array for matrix diagonalization */
	
	cout << "Solving subspace \n";
	
	cout << "Subspace is: " << sub_space_size << "\n";
	
	h=dmatrix(1,sub_space_size,1,sub_space_size);
	d=dvector(1,sub_space_size);
	e=dvector(1,sub_space_size);
	
	cout << "h, d and e initialized \n";
	
	for(i=1; i<=sub_space_size; i++)
		for(j=1; j<=sub_space_size; j++)
			h[i][j]=A_proj[i-1][j-1];
			
	cout << "h formed \n";
	
/* Diagonalize the hamiltonian matrix */
	tred2(h,sub_space_size,d,e);
	tqli(d,e,sub_space_size,h);
	
	cout << "Subspace solved \n";
	
/* Sort the eigenvalues in ascending order */	
	for (i=1; i<sub_space_size; i++)
		for (j=i+1; j<=sub_space_size; j++)
			if (d[i]>d[j])
			{
				dummy = d[i];
				d[i] = d[j];
				d[j] = dummy;
				for (k=1; k<=sub_space_size; k++) e[k] = h [k][i];
				for (k=1; k<=sub_space_size; k++) h[k][i] = h[k][j];
				for (k=1; k<=sub_space_size; k++) h[k][j] = e[k];
			}
	// h holds eigen vectors, d holds eigen values
	
	cout << "Eigenvalues sorted \n";
	
	for(i=0; i<params.N_eig; i++)
	{
		lambda[i] = d[i+1];
		for(j=0; j<params.N_eig; j++)
		{
			eig_vec[i][j]=h[i+1][j+1];
			//cout << eig_vec[i][j] << "\t";
		}
		//cout << "\n";
	}
	
	cout << "eig_vec and lambda formed \n";
	
}
    
    void print_current_step(ostream& output)
    {
	int i;
	output << "Current norms of residuals \n";
	for(i=0;i<params.N_eig;i++)
		output << norma(residual[i]) << endl;
	
	output << "Printing guess vectors \n";
	print_matrix(v, output);
	
	output << "Printing sigma vectors \n";
	print_matrix(sigma_v, output);
	
	output << "Printing projected matrix \n";
	print_matrix(A_proj, output);
	
	output << "Printing current eigen vectors \n";
	print_matrix(eig_vec, output);
	
	output << "Printing current eigenvalues \n";
	print_vec(lambda, output);
	
}

    void print_vec(const std::vector<T> &a, ostream& output) 
	{
  		int i;
  
  		for(i=0; i<a.size(); i++)
  			output << a[i] << " ";
  		output << endl;
	}
	
    void print_matrix( const vector<vector<T> > &A, ostream& output) 
	{
  	int i, j;
  
  	for(i=0; i<A.size(); i++)
  	{
  		for(j=0; j<A[i].size(); j++)
  			output << A[i][j] << " ";
  		output << endl;
  	}
}
		
	Matrix<T> Hmltn;
	Davidson_parameters params;
    vector<vector<T> > v;
    vector<vector<T> > sigma_v;
    vector<vector<T> > A_proj;
    vector<T> lambda;
    vector<vector<T> > eig_vec;
    vector<vector<T> > residual;
	vector<T> diag;
    size_t sub_space_size;
    bool converged;
    
};




#endif	/* DAVIDSON_H */