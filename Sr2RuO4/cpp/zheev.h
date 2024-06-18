/*zheev(mat,4,eval,evecs);
wrapper for zheev - diagonalizes Hermitian matrices. adapted from Scott Shaw's zgeev wrapper by RG on Sep 1, 2010. 
--------------
returns eigenvalues (and if asked for, corresponding eigenvectors) in ascending order.
--------------
define matrices as:
complex<double> **h, **evecs;
double *evals;
evals = new double[dim];
h = new complex<double> *[dim];
evecs = new complex<double> *[dim];
for(int i=0;i<dim;i++) {
    h[i]= new complex<double>[dim];
    evecs[i]= new complex<double>[dim];
  }
---------------
if you want eigenvectors, call zheev as 
zheev(h,dim,evals,evecs);
---------------
if not, call as
zheev(h,dim,evals);
----------------
output is such that
evecs^dg.evecs = identity(dim)
evecs^dg.h.evecs gives the eigenvalues in the diagonals, and zeros in the off-diagonals.
(where matrix multiplication is defined as 
adotb[i][j]=0.; for(u=0;u<dim;u++) adotb[i][j]+=a[i][u]*b[u][j];
*/

#include <cmath>
#include <complex>
#include <iostream>
#include <stdlib.h>

using namespace std;
void zheev(complex<double> **H, int n, double *E);
void zheev(complex<double> **H, int n, double *E,
	   complex<double> **Evecs);

complex<double> *zheev_ctof(complex<double> **in, int rows, int cols);
void zheev_ftoc(complex<double> *in, complex<double> **out,
		int rows, int cols);

void zheev_sort(int n, double *E);
void zheev_sort(int n, double *E, complex<double> **Evecs);

extern "C" void zheev_(char *jobz, char *uplo, 
		       int *n, complex<double> *a, int *lda,double* w,
		       complex<double> *work, int *lwork, double *rwork,
		       int *info);

void zheev(complex<double> **H, int n, double *E)
{
  char jobz, UPLO;
  int lda, lwork, info,n_rwork;
  double *rwork;
  complex<double> *a, *work;
  
  jobz = 'N'; /* N/V to not calculate/ calculate eigenvectors
		  of the matrix H.*/
  UPLO = 'U'; /* L/U to store lower/upper half of matrix */

  lda = n; // The leading dimension of the matrix a.
  a = zheev_ctof(H, n, lda); /* Convert the matrix H from double pointer
				C form to single pointer Fortran form. */

  /* Whether we want them or not, we need to define the matrices
     for the eigenvectors, and give their leading dimensions.
     We also create a vector for work space. */
  lwork = 2*n+1;
  work = new complex<double>[lwork];
  n_rwork=3*n-2;
  rwork = new double[n_rwork];

  zheev_( &jobz,&UPLO,  
	  &n, a, &lda, E, 
          work, &lwork, rwork, &info );

  if(info>0) { cout<<"Error: the algorithm failed to converge;"<<info<<" off-diagonal elements of an intermediate tridiagonal form did not converge to zero"<<endl; 
		exit(1);}
  else if(info<0) {cout<<"Error:the "<<info<<"-th argument had an illegal value"<<endl; exit(1);}

  zheev_sort(n,E);

  // Clean up the memory usage

  delete a;
  delete work;
  delete rwork;
}


void zheev(complex<double> **H, int n, double *E,
	   complex<double> **Evecs)
{
  char jobz, UPLO;
  int lda, lwork, info,n_rwork;
  double *rwork;
  complex<double> *a, *work;

  jobz = 'V'; /* N/V to not calculate/ calculate eigenvectors
		  of the matrix H.*/
  UPLO = 'U'; /* L/U to store lower/upper half of matrix */

  lda = n; // The leading dimension of the matrix a.
  a = zheev_ctof(H, n, lda); /* Convert the matrix H from double pointer
				C form to single pointer Fortran form. */
  lwork = 2*n+1;
  work = new complex<double>[lwork];
  n_rwork=3*n-2;
  rwork = new double[n_rwork];
 
  zheev_( &jobz,&UPLO,  
	  &n, a, &lda, E, 
          work, &lwork, rwork, &info );

   if(info>0) { cout<<"Error: the algorithm failed to converge;"<<info<<" off-diagonal elements of an intermediate tridiagonal form did not converge to zero"<<endl;
		exit(1);}
   else if(info<0) {cout<<"Error:the "<<info<<"-th argument had an illegal value"<<endl; exit(1);}

   zheev_ftoc(a, Evecs, lda, n); //check if this is right. perhaps, rows and columns are interchanged in the eigenvectors stored in a. 

  // Now I need to get the eigenvectors into the output array

  // zheev_sort(n, E, Evecs); // Uncomment if soring is needed!!!!!

  // Clean up the memory usage

  delete a;
  delete work;
  delete rwork;
}


complex<double>* zheev_ctof(complex<double> **in, int rows, int cols)
{
  complex<double> *out;
  int i, j;

  out = new complex<double>[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
  return(out);
}


void zheev_ftoc(complex<double> *in, complex<double> **out,
		 int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*cols];
}


// Sort the eigenenergies by their real components.

void zheev_sort(int n, double *E)
{
  double temp;
  double min, v;
  int imin, ct;

  for (int j=0; j<n-1; j++) {
    min = E[j];
    imin = j;
    for (int i=j+1; i<n; i++) {
      if (E[i]<min) {
	min = E[i];
	imin = i;
      }
    }
    if (imin!=j) {
      temp = E[imin];
      E[imin] = E[j];
      E[j] = temp;
    }
  }
}


void zheev_sort(int n, double *E, complex<double> **Evecs)
{//sorting according to eigenvalue  
  double temp;
  complex<double> ctemp;
  double min, v;
  int imin, ct;

  for (int j=0; j<n-1; j++) {
    min = E[j];
    imin = j;
    for (int i=j+1; i<n; i++) {
      if (E[i]<min) {
	min = E[i];
	imin = i;
      }
    }
    if (imin!=j) {
      temp = E[imin];
      E[imin] = E[j];
      E[j] = temp;

      for (int i=0; i<n; i++) {
	ctemp = Evecs[i][imin];
	Evecs[i][imin] = Evecs[i][j];
	Evecs[i][j] = ctemp;
      }
    }
  }
}
