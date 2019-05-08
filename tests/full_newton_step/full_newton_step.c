/* Romeo Orsolino
 *
 * compile using: gcc full_newton_step.c -o full_newton_step -lm */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOLR 1.0e-7
#define DELTA 1.0e-6
#define MAX_STEP 100.0

double compute_diff(double (*f)(double), double y0)
{
	double gradient;
//    double x1 = x0 - delta;
//    double x2 = x0 + delta;
//    double f1 = f(x1,y0);
//    double f2 = f(x2,y0);
//    gradient[0] = (f2 - f1) / (x2 - x1);

    double y1 = y0 - DELTA;
    double y2 = y0 + DELTA;
    double f1 = f(y1);
    double f2 = f(y2);
//	printf("%lf\n", f1);
//	printf("%lf\n", f2);
    gradient = (f2 - f1) / (y2 - y1);
//	printf("%lf\n", gradient);
    return gradient;
}

double fun(double y)
{
	return 0.5*log(y) + y*y;
}

int main ()
{
	double x0, y0, g;
	double grad;
	double der;
	double (*fun_ptr)(double) = &fun;

	x0 = 0.0;
	y0 = -200.0;
	double semi_interval = 1.0e+7;

	double y_max = y0 + semi_interval;
	double y_min = y0 - semi_interval;


	double r = 10000.0;
	double delta_w = 0.0;
	double w = y0;
	while( fabs(r) > TOLR) {
//		printf("%lf\n", w);
		grad = compute_diff((*fun_ptr), w);
		if (fabs(grad)> 1.0e+4) grad = MAX_STEP;
		printf("%lf\n", grad);
//		r = (*fun_ptr)(w);
		r = fun(w);
		delta_w = -r/grad;

//		printf("%lf\n", fabs(r));


//		printf("%lf\n", delta_w);
		w += delta_w;
	}

	double int_x = -w*w;
	double int_y = w;
	printf("the root is:\n");
	printf("%lf\n", int_x);
	printf("%lf\n", int_y);
	return 0;
}

