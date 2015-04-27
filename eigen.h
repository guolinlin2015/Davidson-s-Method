#ifndef EIGEN_H
#define EIGEN_H
/*******************************************************************************
Eigenvalue solvers, tred2 and tqli, from "Numerical Recipes in C" (Cambridge
Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"

#define NR_END 1
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

/******************************************************************************/
void jacobi(float **a, int n,float d[],float **v,int *nrot)
/******************************************************************************
Computes all eigenvalues and eigenvectors of a real symmetric matrix a[1..n][1..n]. 
On output, elements of a above the diagonal are destroyed. 
d[1..n] returns the eigenvalues of a. v[1..n][1..n] is a matrix whose columns contain, 
on output, the normalized eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
*******************************************************************************/
{
	int j,iq,ip,i;
	float tresh, theta, tau, t, sm,s,h,g,c,*b,*z;

	b = vector(1, n);
	z = vector(1, n);
	for(ip=1;ip<=n;ip++){
		for(iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip] = 1.0;
	}
	for (ip=1;ip<=n;ip++){
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}
	*nrot = 0;
	for(i=1;i<=50;i++){
		sm = 0.0;
		for(ip=1;ip<=n-1;ip++){
			for(iq=ip+1;iq<=n;iq++)
				sm += fabs(a[ip][iq]);
		}
		if (sm == 0.0){
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i<4)
			tresh = 0.2*sm/(n*n);
		else
			tresh = 0.0;
		for(ip=1;ip<=n-1;ip++){
			for(iq=ip+1;iq<=n;iq++){
				g=100.0*fabs(a[ip][iq]);
				if(i>4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip]) && (float)(fabs(d[iq])+g)==(float)fabs(d[iq]))
					a[ip][iq] = 0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if((float)(fabs(h)+g) == (float)fabs(h))
						t = (a[ip][iq])/h;
					else{
						theta = 0.5*h/(a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -=h;
					z[iq] +=h;
					d[ip] -=h;
					d[iq] +=h;
					a[ip][iq] = 0.0;
					for(j=1;j<=ip-1;j++){
						ROTATE(a, j, ip, j, iq)
					}
					for(j=ip+1;j<=iq-1;j++){
						ROTATE(a, ip, j, j, iq)
					}
					for(j=iq+1;j<=n;j++){
						ROTATE(a, ip, j, iq, j)
					}
					for(j=1;j<=n;j++){
						ROTATE(v, j, ip, j, iq)
					}
					++(*nrot);
				}
			}
		}
		for(ip=1;ip<=n;ip++){
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip] = 0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}

/******************************************************************************/
void eigsrt(float d[],float **v,int n)
/******************************************************************************
Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi (§11.1) or tqli (§11.3), 
this routine sorts the eigenvalues into descending order, 
and rearranges the columns of v correspondingly. 
The method is straight insertion.
*******************************************************************************/
{
	int k,j,i;
	float p;

	for(i=1;i<n;i++){
		p=d[k=i];
		for(j=i+1;j<=n;j++)
			if (d[j]>=p) p=d[k=j];
		if (k!=i){
			d[k] = d[i];
			d[i] = p;
			for(j=1;j<=n;j++){
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}


/******************************************************************************/
void tred2(float **a, int n, float d[], float e[])
/*******************************************************************************
Householder reduction of a real, symmetric matrix a[1..n][1..n]. 
On output, a is replaced by the orthogonal matrix Q effecting the
transformation. d[1..n] returns the diagonal elements of the tridiagonal matrix,
and e[1..n] the off-diagonal elements, with e[1]=0. Several statements, as noted
in comments, can be omitted if only eigenvalues are to be found, in which case a
contains no useful information on output. Otherwise they are to be included.
*******************************************************************************/
{
	int l,k,j,i;
	float scale,hh,h,g,f;

	for (i=n;i>=2;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 1) {
			for (k=1;k<=l;k++)
				scale += fabs(a[i][k]);
			if (scale == 0.0) /* Skip transformation. */
				e[i]=a[i][l];
			else {
				for (k=1;k<=l;k++) {
					a[i][k] /= scale; /* Use scaled a's for transformation. */
					h += a[i][k]*a[i][k]; /* Form £m in h. */
				}
				f=a[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g; /* Now h is equation (11.2.4). */
				a[i][l]=f-g; /* Store u in the ith row of a. */
				f=0.0;
				for (j=1;j<=l;j++) {
					/* Next statement can be omitted if eigenvectors not wanted */
					a[j][i]=a[i][j]/h; /* Store u/H in ith column of a. */
					g=0.0; /* Form an element of Aáu in g. */
					for (k=1;k<=j;k++)
						g += a[j][k]*a[i][k];
					for (k=j+1;k<=l;k++)
						g += a[k][j]*a[i][k];
					e[j]=g/h; /* Form element of p in temporarily unused element of e. */
					f += e[j]*a[i][j];
				}
				hh=f/(h+h); /* Form K, equation (11.2.11). */
				for (j=1;j<=l;j++) { /* Form q and store in e overwriting p. */
					f=a[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=1;k<=j;k++) /* Reduce a, equation (11.2.13). */
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		} else
			e[i]=a[i][l];
		d[i]=h;
		}
		/* Next statement can be omitted if eigenvectors not wanted */
		d[1]=0.0;
		e[1]=0.0;
		/* Contents of this loop can be omitted if eigenvectors not
		   wanted except for statement d[i]=a[i][i]; */
		for (i=1;i<=n;i++) { /* Begin accumulation of transformation matrices. */
			l=i-1;
		if (d[i]) { /* This block skipped when i=1. */
			for (j=1;j<=l;j++) {
				g=0.0;
				for (k=1;k<=l;k++) /* Use u and u/H stored in a to form PáQ. */
					g += a[i][k]*a[k][j];
				for (k=1;k<=l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i]=a[i][i]; /* This statement remains. */
		a[i][i]=1.0; /* Reset row and column of a to identity matrix for next iteration. */
		for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
	}
}

/******************************************************************************/
void tqli(float d[], float e[], int n, float **z)
/*******************************************************************************
QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors
of a real, symmetric, tridiagonal matrix, or of a real, symmetric matrix
previously reduced by tred2 sec. 11.2. On input, d[1..n] contains the diagonal
elements of the tridiagonal matrix. On output, it returns the eigenvalues. The
vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix, with
e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues,
several lines may be omitted, as noted in the comments. If the eigenvectors of
a tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as the
identity matrix. If the eigenvectors of a matrix that has been reduced by tred2
are required, then z is input as the matrix output by tred2. In either case,
the kth column of z returns the normalized eigenvector corresponding to d[k].
*******************************************************************************/
{
	float pythag(float a, float b);
	int m,l,iter,i,k;
	float s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i]; /* Convenient to renumber the elements of e. */
	e[n]=0.0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l;m<=n-1;m++) { /* Look for a single small subdiagonal element to split the matrix. */
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((float)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 30) printf("Too many iterations in tqli");
				g=(d[l+1]-d[l])/(2.0*e[l]); /* Form shift. */
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); /* This is dm - ks. */
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { /* A plane rotation as in the original QL, followed by Givens */
					f=s*e[i];          /* rotations to restore tridiagonal form.                     */
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) { /* Recover from underflow. */
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					/* Next loop can be omitted if eigenvectors not wanted */
					for (k=1;k<=n;k++) { /* Form eigenvectors. */
						f=z[k][i+1];
						z[k][i+1]=s*z[k][i]+c*f;
						z[k][i]=c*z[k][i]-s*f;
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

/******************************************************************************/
float pythag(float a, float b)
/*******************************************************************************
Computes (a2 + b2)1/2 without destructive underflow or overflow.
*******************************************************************************/
{
	float absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

#endif
