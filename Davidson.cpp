#include <iostream>
#include <ctime>
#include <cstdlib>
#include "QRdecomp.h"
#include "mvutil.h"
#include "eigen.h"
#include "nrutil.h"

int main(int argc,const char* argv[])
{
	long n = 500;
	float tol = 1e-8;
	long k = 4;
	long l = 2*k;
	long m;
	long i,j,ii;
	float sparsity = 0.0000001;

	srand((unsigned)time(NULL));

	//Construct Matrix
	float **A = matrix(1,n,1,n);
	for(i=1;i<=n;i++){
		A[i][i] = i+1;
	}
	for(i=1;i<=n;i++){
		for(j=1;j<=n;j++)
			A[i][j] = A[i][j] + sparsity * rand();
	}
	float **AT = matrix(1,n,1,n);
	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			AT[i][j] = A[j][i];

	for(i=1;i<=n;i++)
		for(j=1;j<=n;j++)
			A[i][j] = (A[i][j]+AT[i][j])/2;

	//Construct Guess Vector
	float **B = eye(1,n,1,l);
	float *c = vector(1,n);

	float eigValue; 

	for (m = l;m <= n;m++){
		//Gram-schedit Orthogonalize
		orthogonal(B,c,1,n,1,m);
		//Subspace matrix
		float **Abar = matrix(1,m,1,m);
		for(i=1;i<=m;i++){
			float *Bi = MatrixGetRow(B, 1, n, i);
			for(j=1;j<=m;j++){
				float *Bj = MatrixGetRow(B, 1, n, j);
				Abar[i][j] = dotProduct(Bi,MVM(A,Bj,1,n,1,m),1,m);
			}
		}

		float *d = vector(1, m);
		float *e = vector(1, m);
		tred2(Abar, m, d, e);
		float **z = matrix(1, m, 1, m);
		tqli(d, e, m, z);

		float eigenValue = d[m];
		float *eigenVector = vector(1, m);
		for(i=1;i<=m;i++)
			eigenVector[i] = z[i][m]; 
		//Form error vector
		float *Qm = vector(1, n);
		for(i=1;i<=m;i++){
			float *Bi = MatrixGetRow(B, 1, n, i);
			float *ABi = MVM(A, Bi, 1, n, 1, n);
			float *aBi = ScaleVector(Bi, eigenValue, 1, n);
			Qm = VectorAdd(Qm, ScaleVector(VectorSubtract(ABi, aBi, 1, n), eigenVector[i], 1, n), 1, n);
		}

		//Weinstein Lower Bound 
		float norm = VectorNorm(Qm, 1, n);
		if (norm < tol){
			eigValue = eigenValue;
			break;
		}

		float *epth = vector(1,n);
		for(i=1;i<=n;i++)
			epth[i] = Qm[i] * (1.0/(eigenValue-A[i][i]));

		float **mat = eye(1,n,1,n);
		for(i=1;i<=m;i++)
			mat = MatrixMatrixM(mat, MatrixSubtract(eye(1, n, 1, n), VectorVectorTM(MatrixGetRow(B, 1, n, i),MatrixGetRow(B, 1, n, i),1,n,1,n),1,n,1,n),1,n,1,n);
		float *bm_new = MVM(mat, epth, 1, n, 1, n);
		VectorNormalize(bm_new, 1, n);

		B = MatrixAddCol(B,bm_new,1,n,1,m);
	}
	std::cout << eigValue << std::endl;
	return 0;
}