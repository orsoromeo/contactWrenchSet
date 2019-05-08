/* Romeo Orsolino
 *
 * With respect to exercise1.c this program tries to solve a system of nonlinear equations using Lapack and BLAS routines.
 * In particular I used the dgemv_ routine to compute the matrix-vector product and I used the dgetrf_ and dgetri_
 * routines to compute the inverse of a 2x2 matrix.
 *
 * Please note that in this file, as well as in exercise1.c and exercise3.c the memory allocation has not been optimized.
 * This means that the code might give origin to memory leaks, but it has been verified as regards segmentation faults.
 *
 * compile using: gcc -Wall exercise2.c -o exercise2 -lm -llapack -lblas */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

void full_newton_step(double* A1, double* w, double* u, int size)
{
/*	This function performs the following steps:
 * 1) transposition of a squared 2x2 matrix and allocation of its values in a 1 dimensional array. This
 * is required by the fact that Fortran uses the column-major order while C uses the row-major order
 * 2) computation of the matrix-vector product required to compute the magnitude the Newton step
 * 3) inversion of the sign of the Newton step, necessary to have the error decrease with the correct behavior.
 */
	double A[size][size];
	A[0][0]= A1[0];  A[0][1]= A1[1]; 	/* matrix A */
	A[1][0]= A1[2];  A[1][1]= A1[3];

	double S[size*size];

	int k,j;
	for (k=0; k<size; k++)		/* to call a Fortran routine from C we */
	{				/* have to transform the matrix */
	  for(j=0; j<size; j++) {
		  S[j+size*k]=A[j][k];
	  }
	}


    char T = 'N'; // y := alpha*A*x + beta*y.
    int m = size;
    int L = size;
    int inc_x = 1.0;
    int inc_y = 1.0;
    double alpha = 1.0;
    double beta = 0.0;
    double *y = malloc((m-1)*sizeof(*y));

    double x[] = {w[0], w[1]};

    dgemv_(&T, &m, &size, &alpha, S, &L, &x, &inc_x, &beta, y, &inc_y); //matrix-vector multiplication with blas
    int i;
    for(i= 0; i<size; ++i){
        u[i]= -y[i];
        }
    free(y);

}

void inverse(double* A, int N)
{
/*	Compute the inverse of a NxN matrix using the corresponding lapack routine */

    double *ipiv = malloc((N+1)*sizeof(*ipiv));

    int lwork = N*N;

    double *work = malloc(lwork*sizeof(*work));

    int INFO;

    dgetrf_(&N,&N,A,&N,ipiv,&INFO);
    dgetri_(&N,A,&N,ipiv,work,&lwork,&INFO);

    free(ipiv);
    free(work);

}

double * differentiate(double (*f)(double, double), double x0, double y0)
{
/* this method performs the differentiation of a function f with respect to all its variables x and y.
 * The differentiation is performed using a numerical difference centered in the point (x0, y0).
 */
	double * difference = malloc(2*sizeof(double));
    const double delta = 1.0e-6; // or similar
    double x1 = x0 - delta;
    double x2 = x0 + delta;
    double f1 = f(x1,y0);
    double f2 = f(x2,y0);
    difference[0] = (f2 - f1) / (x2 - x1);

    double y1 = y0 - delta;
    double y2 = y0 + delta;
    f1 = f(x0,y1);
    f2 = f(x0,y2);
    difference[1] = (f2 - f1) / (x2 - x1);
    return difference;
}

double** createArray(int n, int m)
{
	/*this routine allocates the amount memory needed to the 2d array */
    double* values = calloc(m, sizeof(double));
    double** rows = malloc(n*sizeof(double*));
    int i;
    for (i=0; i<n; ++i)
    {
        rows[i] = values;
    }
    return rows;
}

void destroyArray(double** arr)
{
    free(*arr);
    free(arr);
}

double fun1(double x, double y)
{

//		return y+2.0*x;  // test 1
//		return y-(x-2.0)*(x-2.0);  // test 2
	return exp(2.0*x) - y;  // test 3
}

double fun2(double x, double y)
{
//	return y-x-3;  // test 1
//	return y-(x+2.0)*(x+2.0);  // test 2
	return y*y + x;  // test 3
}

double** change_sign(double** matrix){

	double** M = createArray(2,2);
	double a = matrix[0][0];
	double b = matrix[0][1];
	double c = matrix[1][0];
	double d = matrix[1][1];

	double * row1 = malloc(2*sizeof(double));
	double * row2 = malloc(2*sizeof(double));
	row1[0] = -a; // a
	row1[1] = -b; // b
	row2[0] = -c; // c
	row2[1] = -d; // d
	M[0] = row1;
	M[1] = row2;

	return M;
}

double** compute_gradient(double (*f1)(double, double), double (*f2)(double, double), double x0, double y0){

	double ** gradient = createArray(2,2);
	double * diff1;
	double * diff2;
	double (*fun1_ptr)(double, double) = f1;
	diff1 = differentiate((*fun1_ptr),x0, y0);

	double (*fun2_ptr)(double, double) = f2;
	diff2 = differentiate((*fun2_ptr),x0, y0);

    gradient[0] = diff1;
    gradient[1] = diff2;

	return gradient;
}

int main ()
{
	/* this is the main routine of the code */

	clock_t start = clock();
	const int size = 2;
	double x0, y0;
	double ** grad = createArray(2,2);

	double (*fun1_ptr)(double, double) = &fun1;
	double (*fun2_ptr)(double, double) = &fun2;


	x0 = -0.5;
	y0 = 0.5;
	double r = 1000.0;
	double r_old = r;

	double * w = malloc(2*sizeof(double));
	w[0] = x0; w[1] = y0;
	double * delta_w = malloc(2*sizeof(double));
	double * f = malloc(2*sizeof(double));
	int j = 0;
	double tol = 10e-6;
	int max_iter_num = 1000;

    double A [2*2];

	while( r > tol) {

		grad = compute_gradient((*fun1_ptr),(*fun2_ptr),w[0], w[1]);
		A[0] = grad[0][0]; 		A[1] = grad[0][1];
		A[2] = grad[1][0];		A[3] = grad[1][1];

	    inverse(A, size);

		f[0] = fun1(w[0],w[1]);
		f[1] = fun2(w[0],w[1]);
		r_old = r;
        r = sqrt(f[0]*f[0] + f[1]*f[1]);

        full_newton_step(A,f,delta_w,size);

        if(r > r_old){
        	printf("flipping the direction of the search \n");
        	delta_w[0] = -delta_w[0];
        	delta_w[1] = -delta_w[1];
        }

		w[0]+= delta_w[0];
        w[1]+= delta_w[1];

        j++;

    	if (j>max_iter_num) {
    		printf(" maximum number of iterations exceeded! \n");
    		break;
    	}

	}
	printf("====>result:<====\n");
	printf("x = %lf\n", w[0]);
	printf("y = %lf\n", w[1]);
	printf("final number of iterations:\n");
	printf("j = %i\n", j);
	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf ("seconds: %f  \n", seconds);

	destroyArray(grad);
	return 0;
}

