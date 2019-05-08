/* Romeo Orsolino
 *
 * Please note that this code (as much as exercise1.c and exercise2.c) is finding the solution of the
 * system of nonlinear equations given in the quiz only if the initial guess is near enough to the final solution.
 * This means that it will not find the correct solution for the suggested starting point (0.0, 0.0).
 * A number of possible way could have been tested to make the algorithm converge to the correct solution (like applying only
 * a small percentage of the full Newton step) but I unfortunately had not enough time to perform this.
 *
 * However the code has been written in such a way to adapt it to a system of an arbitrary number of nonlinear equations.
 * This code was designed for a system of N nonlinear equations and N variables.
 *
 * compile using: gcc -Wall exercise3.c -o exercise3 -lm -llapack -lblas */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


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

void full_newton_step(double* A1, double* w, double* u, int size)
{

	double A[size][size];
//	A[0][0]= A1[0];  A[0][1]= A1[1]; 	/* matrix A */
//	A[1][0]= A1[2];  A[1][1]= A1[3];

	{int k,j;
	for (k=0; k<size; k++)		/* to call a Fortran routine from C we */
	{				/* have to transform the matrix */
	  for(j=0; j<size; j++) {
		  A[j][k] = A1[k+size*j];
	  }
	}}

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

    //y = alpha*A*x + beta*y
    dgemv_(&T, &m, &size, &alpha, S, &L, &x, &inc_x, &beta, y, &inc_y); //matrix-vector multiplication with blas
    int i;
    for(i= 0; i<size; ++i){
        u[i]= -y[i];
        }
    free(y);

}

void inverse(double* A, int N)
{
    double *ipiv = malloc((N+1)*sizeof(*ipiv));
    int lwork = N*N;
    double *work = malloc(lwork*sizeof(*work));
    int INFO;
    dgetrf_(&N,&N,A,&N,ipiv,&INFO);
    dgetri_(&N,A,&N,ipiv,work,&lwork,&INFO);
    free(ipiv);
    free(work);
}

void differentiate(double * (*f)(double *, const int),const int size, double gradient[size][size], double* w0)
{

//	double gradient[size][size];// = malloc(size*size*sizeof(double));
    const double delta = 1.0e-6; // or similar
    double x[size];// = malloc(size*sizeof(double));
    double x_delta1, x_delta2;
    int i,k;


    for(k=0; k<size; k++){

    	for(i=0; i<size; i++){
    		x[i] = w0[i];
    	}

    	for(i=0; i<size; i++){
    		x_delta1 = w0[i] - delta;
    		x_delta2 = w0[i] + delta;
    		x[i] = w0[i] - delta; // x[i] = (*(x + i)) ---- x = &x[0]
    	    double * f1;// = malloc(size*sizeof(double));
    		f1 = f(x,size);
    		double fvalue1 = f1[k];
    		x[i] = w0[i] + delta;
    		double * f2;// = malloc(size*sizeof(double));
    		f2 = f(x,size);
    		double fvalue2 = f2[k];
    		double tmp = (fvalue2 - fvalue1) / (x_delta2 - x_delta1);

    		gradient[k][i] = tmp;
    	}

    }

}


void destroyArray(double** arr)
{
    free(*arr);
    free(arr);
}

double * fun(double * x, const int size)
{
	/* here is where the user should write the nonlinear equations that belong to the system.
	 * Please note that the array func should have the same number of elements that is indicated by the parameter "size"
	 * The user should therefore change the value of "size" in the main function to the desired number of
	 * nonlinear equation that belong to the system.
	 *
	 * Below you can find 3 possible systems of equations that can be tested separately */
	double * func = malloc(size*sizeof(double));

	//////////// Case with 1 equation  //////////////////
	/* Test 1 (linear problem) */
//	func[0] = x[1]+2.0*x[0];

	//////////// Case with 2 equations  //////////////////
	/* Test 1 (linear problem) */
//	func[0] = x[1]+2.0*x[0];
//	func[1] = x[1]-x[0]-3.0;

	/* Test 2  (convex problem) */
//	func[0] = x[1]-(x[0]-2.0)*(x[0]-2.0);
//	func[1] = x[1]-(x[0]+2.0)*(x[0]+2.0);

	/* Test 3 (non-convex problem) */
	func[0] = exp(2.0*x[0]) - x[1];
	func[1] = x[1]*x[1] + x[0];

	//////////// Case with 3 equations  //////////////////
	/* Test 1 (linear problem) */
//	func[0] = x[1]+2.0*x[0]+x[2];
//	func[1] = x[1]-x[0]-x[2];
//	func[2] = x[1]-x[0]-3.0-x[2];

	return func;
}

void full_newton_step_method(const int size, double * w0){

	double grad[size][size];// = createArray(size,size);

	double * (*fun_ptr)(double *, const int) = &fun;

	double r = 1000.0;

	double r_old = r;

	double * w = malloc(size*sizeof(double));

	w = w0;

	double * delta_w = malloc(size*sizeof(double));
	double * f = malloc(size*sizeof(double));
	int j = 0;
	double tol = 10e-6;
	int max_iter_num = 1000;

    double A [size*size];

	while( r > tol) {

		differentiate((*fun_ptr), size, grad, w);

		int k, z;
		for(k=0; k<size; k++){
			for(z=0; z<size; z++){
					A[z+size*k] = grad[k][z];
				}
		}
		/* perform the inversion of the matrix with BLAS routine */
	    inverse(A, size);
	    f = fun(w, size);
		r_old = r;
        int i;
        double r_squared = 0.0;
        for(i = 0;i<size; i++){
        	r_squared += f[i]*f[i];
        }
        r = sqrt(r_squared);
//        r = sqrt(f[0]*f[0] + f[1]*f[1]);

        full_newton_step(A,f,delta_w,size);

        if(r > r_old){
        	printf("flipping the direction of the search \n");

            int i;
            for(i = 0;i<size; i++){
            	delta_w[i] = -delta_w[i];
            }
        }

        for(i = 0;i<size; i++){
        	w[i]+= delta_w[i];
        }

        j++;
    	if (j>max_iter_num) {
    		printf(" maximum number of iterations exceeded! \n");
    		break;
    	}

	}
	printf("====>result:<====\n");
    int i;
    for(i = 0;i<size; i++){
    	printf("%lf\n", w[i]);
    }
	printf("final number of iterations:\n");
	printf("j = %i\n", j);

	free(w);
	free(delta_w);
	free(f);

}

int main ()
{
	clock_t start = clock();

	/* the user should set here the number of equations of the system:*/
	const int size = 2;

	double * w0 = malloc(size*sizeof(double));
	/* The user should add here the initial guess w0. Otherwise another possibility
	 * is to set the initial guess equal to the central point of the feasible domain of the state w */
	w0[0] = -0.5;
	w0[1] = 0.5;
	w0[1] = 0.5;

	/* Perform the Full Newton Step Algorithm */
	full_newton_step_method(size, w0);

	clock_t end = clock();
	float seconds = (float)(end - start) / CLOCKS_PER_SEC;
	printf ("seconds: %f  \n", seconds);

	return 0;
}

