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

double besselj_function(double x, void * params);
double besseli_function(double x, void * params);
double inverse_laplace_transform(double t, double k, double alpha, double ksi, double c, double freq, bool flag);


// }