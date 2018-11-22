#include "bessel.h"
#include <iostream>
#include <cmath>


// extern "C" {

double inverse_laplace_transform_hx(double t, double k, double alpha, double ksi, double c, double freq, bool flag)
{
   if (t<=k)
      return 0;
   else
   {
      // numeric integration
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      double result, error;
      //~ size_t nevals;
      gsl_function F;

      // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
      // double c_o = 0.3*c;
      
      struct my_data data;
      data.t = t;
      data.k = k;
      data.alpha = alpha;
      data.ksi = ksi;
      //~ data.c = c;
      data.freq = freq;
      
      if (flag)
         F.function = &besselj_function_hx;
      else
         F.function = &besseli_function_hx;
      
      F.params = (void *)(&data);
      
      // printf ("result = % .18f\n", result);
      // printf ("estimated error = % .18f\n", error);
      // //printf ("intervals = %zu\n", w->size);
      // printf ("Size of interval = %g\n", t-k);
      gsl_integration_qags(&F, k, t, 1e-12, 1e-8, 1000, w, &result, &error);
      //printf ("intervals = %zu\n", w->size);
      gsl_integration_workspace_free(w);
      
      // std::cout << "first term       = " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;         
      // std::cout << "second term        = " << alpha*k*result << std::endl;
      // std::cout << "integral error   = " << error << std::endl;
      
      if (flag)
      {
         return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))-alpha*k*result;
      }
      else
      {
         return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))+alpha*k*result;
      }
   }
}

double besselj_function_hx(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   
   gsl_sf_result res;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   /*int status = */gsl_sf_bessel_J1_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}

double besseli_function_hx(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   gsl_sf_result res;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   /*int status = */gsl_sf_bessel_I1_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}
   
double inverse_laplace_transform_hz(double t, double k, double alpha, double ksi, double c, double freq, bool flag)
{
   if (t<=k)
      return 0;
   else
   {
      // numeric integration
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      double result, error;
      //~ size_t nevals;
      gsl_function F;

      // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
      // double c_o = 0.3*c;
      
      struct my_data data;
      data.t = t;
      data.k = k;
      data.alpha = alpha;
      data.ksi = ksi;
      //~ data.c = c;
      data.freq = freq;
      
      if (flag)
         F.function = &besselj_function_hz;
      else
         F.function = &besseli_function_hz;
      
      F.params = (void *)(&data);
      
      // printf ("result = % .18f\n", result);
      // printf ("estimated error = % .18f\n", error);
      // //printf ("intervals = %zu\n", w->size);
      // printf ("Size of interval = %g\n", t-k);
      gsl_integration_qags(&F, k, t, 1e-12, 1e-8, 1000, w, &result, &error);
      //printf ("intervals = %zu\n", w->size);
      gsl_integration_workspace_free(w);
      
      // std::cout << "first term       = " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;         
      // std::cout << "second term        = " << alpha*k*result << std::endl;
      // std::cout << "integral error   = " << error << std::endl;
      
      return result;
   }
}

double besselj_function_hz(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   
   gsl_sf_result res;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   /*int status = */gsl_sf_bessel_J0_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}

double besseli_function_hz(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   gsl_sf_result res;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   /*int status = */gsl_sf_bessel_I0_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}   
   
double inverse_laplace_transform_ey(double t, double k, double alpha, double ksi, double c, double freq, bool flag)
{
   if (t<=k)
      return 0;
   else
   {
      // numeric integration
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      double result, error;
      //~ size_t nevals;
      gsl_function F;

      // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
      // double c_o = 0.3*c;
      
      struct my_data data;
      data.t = t;
      data.k = k;
      data.alpha = alpha;
      data.ksi = ksi;
      //~ data.c = c;
      data.freq = freq;
      
      if (flag)
         F.function = &besselj_function_ey;
      else
         F.function = &besseli_function_ey;
      
      F.params = (void *)(&data);
      
      // printf ("result = % .18f\n", result);
      // printf ("estimated error = % .18f\n", error);
      // //printf ("intervals = %zu\n", w->size);
      // printf ("Size of interval = %g\n", t-k);
      gsl_integration_qags(&F, k, t, 1e-12, 1e-8, 1000, w, &result, &error);
      //printf ("intervals = %zu\n", w->size);
      gsl_integration_workspace_free(w);
      
      // std::cout << "first term       = " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;         
      // std::cout << "second term        = " << alpha*k*result << std::endl;
      // std::cout << "integral error   = " << error << std::endl;
      
      // if (flag)
      // {
         // return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))-alpha*result;
      // }
      // else
      // {
         // return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))+alpha*result;
      // }
      
      return result;
   }
}

double besselj_function_ey(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   // std::cout << "miao" << std::endl;
   gsl_sf_result res;//, res2;
   // std::cout << "miao" << (x/k) << std::endl;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   
   /*int status = */gsl_sf_bessel_J0_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*(2*PIE*freq)*cos(2*PIE*freq*(t-x));
   
   // int status   = gsl_sf_bessel_J1_e(argument,&res)/sqrt( pow(x,2) - pow(k,2) );
   // double f     = exp(-ksi*x)*x*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));

   
   return f;
}

double besseli_function_ey(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // std::cout << "ciao" << std::endl;
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   gsl_sf_result res;//, res2;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   

   /*int status = */gsl_sf_bessel_I0_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*2*PIE*freq*cos(2*PIE*freq*(t-x));
   
   
   // int status       =    gsl_sf_bessel_I1_e(argument,&res)/sqrt( pow(x,2) - pow(k,2) );
   // double f        = exp(-ksi*x)*x*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}   

double inverse_laplace_transform_ey_old(double t, double k, double alpha, double ksi, double c, double freq, bool flag)
{
   if (t<=k)
      return 0;
   else
   {
      // numeric integration
      gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
      double result, error;
      //~ size_t nevals;
      gsl_function F;

      // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
      // double c_o = 0.3*c;
      
      struct my_data data;
      data.t = t;
      data.k = k;
      data.alpha = alpha;
      data.ksi = ksi;
      //~ data.c = c;
      data.freq = freq;
      
      if (flag)
         F.function = &besselj_function_ey_old;
      else
         F.function = &besseli_function_ey_old;
      
      F.params = (void *)(&data);
      
      // printf ("result = % .18f\n", result);
      // printf ("estimated error = % .18f\n", error);
      // //printf ("intervals = %zu\n", w->size);
      // printf ("Size of interval = %g\n", t-k);
      gsl_integration_qags(&F, k, t, 1e-12, 1e-8, 1000, w, &result, &error);
      //printf ("intervals = %zu\n", w->size);
      gsl_integration_workspace_free(w);
      
      // std::cout << "first term       = " << exp(-ksi*(k))*sin(2*PIE*freq*(t-k)) << std::endl;         
      // std::cout << "second term        = " << alpha*k*result << std::endl;
      // std::cout << "integral error   = " << error << std::endl;
      
      if (flag)
      {
         
         return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))-alpha*k*result;
      }
      else
      {
         return exp(-ksi*(k))*sin(2*PIE*freq*(t-k))+alpha*k*result;
      }
   }
}

double besselj_function_ey_old(double x, void * params)
{
   
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   
   gsl_sf_result res;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   /*int status = */gsl_sf_bessel_J1_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}

double besseli_function_ey_old(double x, void * params)
{
   struct my_data *pt = (my_data *)(params);

   double t    = pt->t;
   double k     = pt->k;
   double alpha = pt->alpha;
   double ksi   = pt->ksi;
   //~ double c     = pt->c;
   double freq  = pt->freq;
   
   // double c_o = 1/sqrt(3*8.854187817e-12*4e-7*3.141592);
   // double c_o   = 0.3*c;
   gsl_sf_result res;
   double argument = alpha*sqrt( pow(x,2) - pow(k,2) );
   /*int status = */gsl_sf_bessel_I1_e(argument,&res);
   double f     = exp(-ksi*x)*res.val*sin(2*PIE*freq*(t-x))/sqrt( pow(x,2) - pow(k,2));

   // double k     = pt->k;
   // double alpha = pt->alpha;
   // double f     = gsl_sf_bessel_J1(alpha*sqrt( pow(x,2) - pow(k,2) ))/ sqrt( pow(x,2) - pow(k,2) );
   
   return f;
}   
// }

