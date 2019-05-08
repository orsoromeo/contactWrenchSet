/* Romeo Orsolino
 *
 * This function computes the root of two nonlinear equations. It has been tested
 * for three different functions (that can be used inside the "fun1" and "fun2" methods by decommenting the
 * desired line).
 *
 * In the case of the system of equations given in the Quiz the root is unfortunately only found
 * when the initial guess is given in a neighbourhood of the final solution. The only attempt I did to
 * make the solution more "global" was by implementing a simple "backtracking" (by simply changing the sign
 * of the Newton direction when this does not decrease the error). Many more refined techniques require could have
 * been implemented: by instance, the use of a parameter alpha whose value depends on how much the error
 * decreased in the previous iteration.
 *
 * compile using: gcc exercise1.c -o exercise1 -lm */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double * product(double** grad, double* r)
{

	double * prod = malloc(2*sizeof(double));
	double a = grad[0][0]; double b = grad[0][1];
	double c = grad[1][0]; double d = grad[1][1];
	double r1 = r[0]; double r2 = r[1];
    prod[0] = a*r1 + b*r2;
    prod[1] = c*r1 + d*r2;

    return prod;
}

double * differentiate(double (*f)(double, double), double x0, double y0)
{

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

	double** transposed = createArray(2,2);
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
	transposed[0] = row1;
	transposed[1] = row2;

	return transposed;
}

double** inverse2d(double** matrix, int m, int n){

	double** inverse = createArray(2,2);
	double a = matrix[0][0];
	double b = matrix[0][1];
	double c = matrix[1][0];
	double d = matrix[1][1];

	double det = a*b - c*d;

	double * row1 = malloc(2*sizeof(double));
	double * row2 = malloc(2*sizeof(double));
	row1[0] = d / det; // a
	row1[1] = - b / det;// b
	row2[0] = -c / det; // c
	row2[1] = a / det;// d
	inverse[0] = row1;
	inverse[1] = row2;

	return inverse;
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
	clock_t start = clock();
	double x0, y0;
	double ** grad = createArray(2,2);

	double** inv = createArray(2,2);
	double ** minus_inv = createArray(2,2);

	double (*fun1_ptr)(double, double) = &fun1;
	double (*fun2_ptr)(double, double) = &fun2;

	x0 = -0.5;
	y0 = 0.5;
	double r = 1000.0;
	double r_old = r;
//	double delta_w = 0.0;
	double * w = malloc(2*sizeof(double));
	w[0] = x0; w[1] = y0;
	double * delta_w = malloc(2*sizeof(double));
	double * f = malloc(2*sizeof(double));
	int j = 0;
	double tol = 10e-6;
	double lambda = 1.0;
	int max_iter_num = 10000;
	while( r > tol) {

		grad = compute_gradient((*fun1_ptr),(*fun2_ptr),w[0], w[1]);

		inv = inverse2d(grad, 2, 2);

		minus_inv = change_sign(inv);

		f[0] = fun1(w[0],w[1]);
		f[1] = fun2(w[0],w[1]);
		r_old = r;
        r = sqrt(f[0]*f[0] + f[1]*f[1]);
		delta_w = product(minus_inv,f);
        if(r > r_old){
        	//	implement a proper backtracking algorithm if you have time
        	printf("Backtracking along the Newton direction \n");
        	lambda = -1.0;
        }

//    	printf("error:\n");
//    	printf("%lf\n", r);

		w[0]+= delta_w[0]*lambda;
        w[1]+= delta_w[1]*lambda;

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

	free(w);
	free(delta_w);
	free(f);
	return 0;
}

