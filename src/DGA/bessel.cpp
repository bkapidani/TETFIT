#include "bessel.h"
#include <iostream>
#include <cmath>


// extern "C" {

	double inverse_laplace_transform(double t, double k, double alpha, double ksi, double c, double freq, bool flag)
	{
		if (t<=k)
			return 0;
		else
		{
			// numeric integration
			gsl_integration_workspace * w = gsl_integration_workspace_alloc(2000);
			double result, error;
			gsl_function F;

			// double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
			// double c_o = 0.3*c;
			
			struct my_data data;
			data.t = t;
			data.k = k;
			data.alpha = alpha;
			data.ksi = ksi;
			data.c = c;
			data.freq = freq;
			
			if (flag)
				F.function = &besselj_function;
			else
				F.function = &besseli_function;
			
			F.params = (void *)(&data);
			
			// printf ("result = % .18f\n", result);
			// printf ("estimated error = % .18f\n", error);
			// printf ("intervals = %zu\n", w->size);
			// printf ("Size of interval = %g\n", t-k);
			gsl_integration_qags(&F, k, t, 0, 1e-8, 2000, w, &result, &error);
			gsl_integration_workspace_free(w);
			
			// std::cout << "first term 		= " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;			
			// std::cout << "second term     	= " << alpha*k*result << std::endl;
			// std::cout << "integral error	= " << error << std::endl;
			
			if (flag)
			{
				// std::cout << "time = " << t << " a = " << k << " b = " << t << std::endl;
				// std::cout << "ksi = " << ksi << " alpha = " << alpha << " value = " << alpha*k*result << std::endl;
				// std::cout << "1st term = " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;
				return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))-alpha*k*result;
			}
			else
			{
				// std::cout << "besseli" << std::endl;
				// if (t>5.9e-10 && t< 6.1e-10)
				// {
					// std::cout << "time = " << t << " a = " << k << " b = " << t << std::endl;
					// std::cout << "ksi = " << ksi << " alpha = " << alpha << " value = " << alpha*k*result << std::endl;
					// std::cout << "1st term = " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;
				// }
				return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))+alpha*k*result;
			}
		}
	}
	
	
	double besselj_function(double x, void * params)
	{
		
        struct my_data *pt = (my_data *)(params);
	
		double t	 = pt->t;
		double k     = pt->k;
		double alpha = pt->alpha;
		double ksi   = pt->ksi;
		double c     = pt->c;
		double freq  = pt->freq;
		
		// double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
		// double c_o   = 0.3*c;
		
		gsl_sf_result res;
		double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
		int status = gsl_sf_bessel_J1_e(argument,&res);
		double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));
 
		// double k     = pt->k;
		// double alpha = pt->alpha;
		// double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
		
		return f;
	}

	double besseli_function(double x, void * params)
	{
		
        struct my_data *pt = (my_data *)(params);
	
		double t	 = pt->t;
		double k     = pt->k;
		double alpha = pt->alpha;
		double ksi   = pt->ksi;
		double c     = pt->c;
		double freq  = pt->freq;
		
		// double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
		// double c_o   = 0.3*c;
		gsl_sf_result res;
		double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
		int status = gsl_sf_bessel_I1_e(argument,&res);
		double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));
 
		// double k     = pt->k;
		// double alpha = pt->alpha;
		// double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
		
		return f;
	}
// }


