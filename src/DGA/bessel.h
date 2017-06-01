/* file bessel.h */
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>
#include <stdlib.h>

const long double PIE = 3.141592653589793238L;

struct my_data
{
	double t, ksi, c, freq;
    double k, alpha;
};

// extern "C" {

double besselj_function_hx(double x, void * params);
double besseli_function_hx(double x, void * params);
double besselj_function_hz(double x, void * params);
double besseli_function_hz(double x, void * params);
double besselj_function_ey(double x, void * params);
double besseli_function_ey(double x, void * params);
double besselj_function_ey_old(double x, void * params);
double besseli_function_ey_old(double x, void * params);
double inverse_laplace_transform_hx(double t, double k, double alpha, double ksi, double c, double freq, bool flag);
double inverse_laplace_transform_hz(double t, double k, double alpha, double ksi, double c, double freq, bool flag);
double inverse_laplace_transform_ey(double t, double k, double alpha, double ksi, double c, double freq, bool flag);
double inverse_laplace_transform_ey_old(double t, double k, double alpha, double ksi, double c, double freq, bool flag);
// }