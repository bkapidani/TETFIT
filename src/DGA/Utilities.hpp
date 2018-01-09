/*
 * This source file is part of FDTD UNIUD.
 *
 * Copyright (C) 2017, Bernard Kapidani - kapidani.bernard@spes.uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
 
 //file Utilities.hpp
#pragma once

#include <iostream>
#include <array>
#include <stdexcept>
#include <stdio.h>
#include <cstdlib>
#include <utility> 
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <map>
#include <iomanip>
#include <set>
#include <thread>
// #include <limits>
#include <sstream>
#include <silo.h>
// #include <future>
#include <chrono>
#include <mutex>

// #ifndef EIGEN_USE_BLAS 
// #define EIGEN_USE_BLAS
// #endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Eigenvalues>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/GenEigsSolver.h> 
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseCholesky.h>


#include "timecounter.h"
#include "sgnint32_t.hpp"
#include "mapped_file.h"
#include "strtot.hpp"
#include "bessel.h"
#include "burkardt.h"
#include "Seb.h" // from Fischer
#include "date.h" //from howard hinnant
#include "agmg.hpp"
#include "Efield.hpp"
#include "Hfield.hpp"

#pragma once

/*general parameters*/
const long double PI = 3.141592653589793238L;
double mu0 = 4*PI*1e-7;
double epsilon0 = 8.854187817e-12;
double c0 = 1 / sqrt( mu0 * epsilon0 );

namespace parser 
{
	template<typename T>
	std::tuple<T, T, T>
	read_point_line(const char *str, char **endptr)
	{
		T t1, t2, t3;
		
		t1 = strtot<T>(str, endptr);
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		
		//return std::make_tuple(t1/1000.0, t2/1000.0, t3/1000.0);
		return std::make_tuple(t1, t2, t3);
	}

	template<typename T>
	std::tuple<T, T, T>
	read_gmsh_point_line(const char *str, char **endptr)
	{
		T t0, t1, t2, t3;
		
		t0 = strtot<T>(str, endptr); 
		t1 = strtot<T>(*endptr, endptr);
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		
		//return std::make_tuple(t1/1000.0, t2/1000.0, t3/1000.0);
		return std::make_tuple(t1, t2, t3);
	}

	template<typename T>
	std::tuple<T, T, T, T, T>
	read_tetrahedron_line(const char *str, char **endptr)
	{
		T t1, t2, t3, t4, t5;
		
		t1 = strtot<T>(str, endptr);
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		t4 = strtot<T>(*endptr, endptr);
		t5 = strtot<T>(*endptr, endptr);
		
		return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
	}
		
	template<typename T>
	std::tuple<T, T, T, T>
	read_triangle_line(const char *str, char **endptr)
	{
		T t1, t2, t3, t4;
		
		t1 = strtot<T>(str, endptr);
		t2 = strtot<T>(*endptr, endptr);
		t3 = strtot<T>(*endptr, endptr);
		t4 = strtot<T>(*endptr, endptr);
		
		return std::make_tuple(t1, t2-1, t3-1, t4-1);
	}

	template<typename T>
	std::tuple<T, T, T, T, T, T>
	read_element_line(const char *str, char **endptr)
	{
		T t0, t1, t2, t3, t4, t5;
		
		// std::ofstream os;
		// os.open("./debug/parser_tets.dat", std::ios::out | std::ios::app );
		
		t0 = strtot<T>(str, endptr);
		t1 = strtot<T>(*endptr, endptr);
		
		if (t1 == 4) //tetrehedron
		{
			t0 = t1;
			t1 = strtot<T>(*endptr, endptr);
			t1 = strtot<T>(*endptr, endptr);
			
			
			t2 = strtot<T>(*endptr, endptr);
			t2 = strtot<T>(*endptr, endptr);
			t3 = strtot<T>(*endptr, endptr);
			t4 = strtot<T>(*endptr, endptr);
			t5 = strtot<T>(*endptr, endptr);
			
			
			// os << t0 << " " << t1 << " " << t2 << " " << t3 << " " << t4 << " " << t5 << std::endl;  
		}
		else if (t1 == 2) //triangle
		{
			t0 = t1;
			t1 = strtot<T>(*endptr, endptr);
			t1 = strtot<T>(*endptr, endptr);
			
			
			t2 = strtot<T>(*endptr, endptr);
			t2 = strtot<T>(*endptr, endptr);
			t3 = strtot<T>(*endptr, endptr);
			t4 = strtot<T>(*endptr, endptr);
			t5 = 1;
		}
		else if (t1 == 1) //physical line
		{
			t0 = t1;
			T n_of_tags = strtot<T>(*endptr, endptr);
			for (size_t k = 0; k < n_of_tags; k++)
			{
				T t_token = strtot<T>(*endptr, endptr);
				if (k==0)
				{
					t1 = t_token;
				}
			}
			
			t2 = strtot<T>(*endptr, endptr);
			t3 = strtot<T>(*endptr, endptr);
			t4 = 1;
			t5 = 1;
		}	
		else
		{
			t0 = t1 = 0;
			t2 = t3 = t4 = t5 = 1;
		}
		
		// os.close();
		return std::make_tuple(t0, t1, t2-1, t3-1, t4-1, t5-1);
	}
	
} //namespace parser


typedef std::vector<sgnint32_t<int32_t>>						cluster_list;
typedef std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>			volume_type;
typedef std::tuple<uint32_t,uint32_t,uint32_t>					surface_type;
typedef std::pair<surface_type,uint32_t>						label_surface_type;
typedef std::tuple<uint32_t,uint32_t>							edge_type;
typedef std::pair<edge_type,uint32_t>							label_edge_type;
typedef std::pair<uint32_t,uint32_t>							label_node_type;
typedef std::tuple< volume_type, uint32_t >						tm_tuple;
typedef std::tuple< surface_type, uint32_t >					sm_tuple;
typedef Eigen::Triplet<double,uint32_t>							double_triplet;
typedef std::string												Primitive;
typedef std::string												Sourcetype;
typedef std::string												Direction;
typedef std::string												BaseFunction;
typedef std::string												BoundaryConditionType;
typedef std::string												Profile;
typedef std::string												Carrier;
typedef std::string												Meshtype;
typedef std::string												OutputMode;
typedef std::string												SimMethod;
typedef std::string												Solver;
typedef std::string												SolidType;
typedef double													Amplitude;
typedef double 													Frequency;
typedef double 													WaveNumber;
typedef double													Duration;
typedef std::array<double,4>									SpaceTimePoint;
typedef std::array<double,3> 									WaveVector;
typedef Seb::Point<double> 										Point;
typedef std::vector<Point> 										PointVector;
typedef Seb::Smallest_enclosing_ball<double> 					Miniball;
		
//lists of allowed string constants
const std::vector<Primitive>							definables		= {"material","source","mesh","bc","simulation","geometry","output","refinement"};
const std::vector<Primitive>							meshdefinables	= {"grid","solid"};
const std::vector<SolidType>							solidtypes		= {"sphere","box","cylinder"};
const std::vector<Sourcetype>   						sourcetypes   	= { "e", "b", "j", "h", "d" };
const std::vector<Profile>   							profiles   		= { "wave", "gaussian", "const" };
const std::vector<Carrier>   							carriers   		= { "sin", "cos", "gaussian", "dc", "ricker" };
const std::vector<Direction>    						directions    	= { "x", "y", "z", "r", "psi", "rho", "phi", "theta" };
const std::vector<Direction>    						bc_directions   = { "+x", "+y", "+z", "r", "+t", "rho", "-t", "-x", "-y", "-z"};
const std::vector<BaseFunction> 						modes 			= { "sin", "cos" };
const std::vector<BoundaryConditionType>				bctypes			= { "pec", "pmc", "pml" };
const std::vector<Meshtype>								meshtypes		= { "tetrahedral", "cartesian", "none"};
const std::vector<Meshtype>								meshers		    = { "netgen", "gmsh", "none"};
const std::vector<OutputMode>							outputmodes		= { "silo", "probepoint", "maxerror", "l2norm"};
const std::vector<SimMethod>							simmethods		= { "fdtd", "dga", "fem", "dgao2", "frac", "gfdtd", "fraco2", "fdtdo2"};
const std::vector<Solver>								solvers			= { "cg", "agmg"};

const std::runtime_error main_missing_file(std::string("Input file missing! Correct use is: \"tetfit input_file\" "));
const std::runtime_error pml_missing(std::string("Sorry, PML not implemented yet, getting there!"));
const std::runtime_error pmc_missing(std::string("Sorry, PMC not implemented yet, getting there!"));
const std::runtime_error bc_unknown_type(std::string("Unrecognized boundary condition type! Available: pec, pmc, pml"));
const std::runtime_error bc_unknown_parameter(std::string("Unrecognized boundary condition parameter! Available: type"));
const std::runtime_error bc_unknown_direction(std::string("Unrecognized direction! Available: +-x, +-y, +-z, r, +-t, rho"));
const std::runtime_error src_unknown_direction(std::string("Unrecognized direction! Available: x, y, z"));
const std::runtime_error src_unknown_type(std::string("Unrecognized source type!") + 
                                          std::string("Available: e (electric field), b (magnetic field), j (current density)"));
const std::runtime_error src_unknown_profile(std::string("Unrecognized source profile! Available: const, wave, gaussian"));
const std::runtime_error src_unknown_carrier(std::string("Unrecognized source carrier! Available: dc, sin, cos, gaussian, ricker"));
const std::runtime_error coordinates_syntax(std::string("coordinates must be inside braces {...,..,..}"));
const std::runtime_error unbalanced_bracket(std::string("unbalanced bracket"));
const std::runtime_error too_many_coords(std::string("Maximum of three coordinates!"));
const std::runtime_error too_few_coords(std::string("Undefined end to list of coordinates!"));
const std::runtime_error src_unknown_bf(std::string("Unrecognized base function! Available: sin, cos"));
const std::runtime_error src_unknown_parameter(std::string("unrecognized parameter for source!"));
const std::runtime_error material_unknown_parameter(std::string("Unrecognized material parameter! Available: epsilon, mu, sigma, chi"));
const std::runtime_error mesh_unknown_type(std::string("Undefined mesh type! Available: tetrahedral, cartesian"));
const std::runtime_error mesh_inexistent_file(std::string("Mesh input file does not exist or it is not in specified path"));
const std::runtime_error mesh_unknown_mesher(std::string("Undefined mesher! Available: netgen, gmsh"));
const std::runtime_error mesh_unknown_parameter(std::string("Undefined mesh parameter! Available: file, name, type, mesher, scalefactor"));
const std::runtime_error ref_unknown_parameter(std::string("Undefined refinement parameter! Available: x, y, z, vector"));
const std::runtime_error sim_unknown_output(std::string("Undefined output mode type! Available: silo, probepoint, maxerror, l2norm"));
const std::runtime_error sim_unknown_method(std::string("Undefined simulation method! Available: fdtd, dga, fem"));
const std::runtime_error sim_unknown_solver(std::string("Unavailable solver! Available: cg, agmg"));
const std::runtime_error sim_unknown_parameter(std::string("Undefined simulation parameter! Available: source, mesh, duration, output"));
const std::runtime_error out_unknown_parameter(std::string("Undefined output parameter! Available: mode, period, probe, name"));
const std::runtime_error sim_invalid_coord_system(std::string("Invalid system of coordinates! Available: cartesian, cylindrical, spherical"));
const std::runtime_error set_wo_define(std::string("define something before setting variables"));
const std::runtime_error unknown_define(std::string("can only define one of the following: material, mesh, refinement, boundary condition, simulation, output, source"));
const std::runtime_error grid_unknown_define(std::string("Unknown primivite: can only define solid inside cartesian meshes"));
const std::runtime_error end_wo_define(std::string("ending unopened definition"));
const std::runtime_error unknown_instruction(std::string("Unknown instruction inside define block"));
const std::runtime_error unexpected_end(std::string("File ended unexpectedly before ending definition"));
const std::runtime_error out_of_bounds_freq(std::string("Output frequency must be a value between 0 and simulation time!"));
const std::runtime_error solid_unknown_parameter(std::string("Unrecognized solid parameter! Available: type, radius, center, corner, size"));
const std::runtime_error solid_unknown_type(std::string("Unavailable solid type! Available: sphere, box, cylinder"));
const std::runtime_error solid_negative_value(std::string("Positive definite quantity forced to negative value!"));
const std::runtime_error incompatible_meth_mesh(std::string("If mesh type is tetrahedral, method must be fem or dga. If mesh type is cartesian, method must be fdtd"));

std::string mesh_throw_preamble("In cartesian mesh definition file: ");

template<typename T>
void sort_unique(std::vector<T>& v) //useful as stand-alone
{
	std::sort(v.begin(), v.end());
	auto uniq_iter = std::unique(v.begin(), v.end());
	v.erase(uniq_iter, v.end());
}

typedef std::string StringType;
StringType GetBaseFilename(const char *filename)
{
    StringType fName(filename);
    size_t pos = fName.rfind(".");
    if(pos == StringType::npos)  //No extension.
        return fName;

    if(pos == 0)    //. is at the front. Not an extension.
        return fName;

    return fName.substr(0, pos);
}

std::string tail(std::string const& source, size_t const length) 
{
	if (length >= source.size())
		return source;
	return source.substr(source.size() - length);
} // tail

class add_to_sparse 
{
	public:
	// add_to_sparse(const double& a, const double& b) : a(a), b(b) {}
	add_to_sparse() {}
	double operator()(const double& a, const double& b) const { return a+b; }
	private:
	// double a, b;
};

class overwrite_to_sparse
{
	public:
	// overwrite_to_sparse(const double& a, const double& b) : a(a), b(b) {}
	overwrite_to_sparse() {}
	double operator()(const double& a, const double& b) const { return b; }
	
	private:
	double a, b;
};

void MyThrow(uint32_t input_line, const std::runtime_error& e)
{
	std::cout << "Input file error at line " << input_line << ": " << e.what() << std::endl;
	throw e;
}

void DateAndTime(void)
{
	using namespace date;
	using namespace std::chrono;
	std::cout << "On " << system_clock::now() << '\n';
}

std::pair<Eigen::Vector3d,Eigen::Vector3d> analytic_value_excite_h(SpaceTimePoint p, double sigma, double eps, double mu, double freq)
{
	auto x = p[0]; auto y = p[1]; auto z = p[2]; auto t = p[3];
	double c = 1/sqrt(eps*mu);
	// double ax=5e-2;
	// double az=10e-2;
	// double ay=2.5e-2; //momentarily useless
	double ax=1;
	double az=1;
	double ay=1; //momentarily useless
	double ksi = sigma/eps/2;
	double alph1 = pow(PI*c/ax,2);
	double alph2 = 0.25*pow(sigma/eps,2);
	double alpha;
	bool flag = true;
	if (alph1>alph2)
	{
		flag = true;
		alpha =  sqrt(alph1-alph2);
	}
	else
	{
		flag = false;
		alpha = sqrt(alph2-alph1);
	}
	
	double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12;
	int32_t j;
	// i=j=0;
	a1=a2=a3=a4=a5=a6=a7=a8=a9=a10=a11=a12=0;
	
	/******************************HX******************************************************/
	j=0;
	while (true)
	{
		double k = (z + 4 * az * j) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			// std::cout << "First set of waves: integrate from " << k << " to " << t << " seconds" << std::endl;
			// std::cout << "k = " << k <<  " time = " << t << std::endl;
			a1 += inverse_laplace_transform_hx(t,k, alpha, ksi, c, freq, flag);
			a5 += inverse_laplace_transform_hz(t,k, alpha, ksi, c, freq, flag);
			// a9 -= inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
			a9 -= inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	j=0;
	while (true)
	{
		double k = (4*az*j + 2*az + z) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			// std::cout << "Second set of waves: integrate from " << k << " to " << t << " seconds" << std::endl;
			// std::cout << "k = " << k <<  " time = " << t << std::endl;
			 a2 -= inverse_laplace_transform_hx(t,k, alpha, ksi, c, freq, flag);
			 a6 -= inverse_laplace_transform_hz(t,k, alpha, ksi, c, freq, flag);
			// a10 += inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
			a10 += inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	j=0;
	while (true)
	{
		double k = (4*az*j + 2*az - z) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			// std::cout << "First set of waves: integrate from " << k << " to " << t << " seconds" << std::endl;
			// std::cout << "k = " << k <<  " time = " << t << std::endl;
			 a3 += inverse_laplace_transform_hx(t,k, alpha, ksi, c, freq, flag);
			 a7 -= inverse_laplace_transform_hz(t,k, alpha, ksi, c, freq, flag);
			// a11 += inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
			a11 += inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	j=0;
	while (true)
	{
		double k = (4*az*j + 4*az - z) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			// std::cout << "Second set of waves: integrate from " << k << " to " << t << " seconds" << std::endl;
			// std::cout << "k = " << k <<  " time = " << t << std::endl;
			 a4 -= inverse_laplace_transform_hx(t,k, alpha, ksi, c, freq, flag);
			 a8 += inverse_laplace_transform_hz(t,k, alpha, ksi, c, freq, flag);
			// a12 -= inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
			a12 -= inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	double hx_double = sin(PI*x/ax)*(a1+a2+a3+a4);
	double hy_double = 0;
	double hz_double = (PI/ax)*c*cos(PI*x/ax)*(a5+a6+a7+a8);
	double ex_double = 0;
	// double ey_double = sqrt(mu/eps)*sin(PI*x/ax)*((a9+a10+a11+a12)-sigma*(a5+a6+a7+a8));
	double ey_double = sqrt(mu/eps)*sin(PI*x/ax)*(a9+a10+a11+a12);
	double ez_double = 0;
	
	return std::make_pair<Eigen::Vector3d,Eigen::Vector3d>(Eigen::Vector3d({ex_double,ey_double,ez_double}),Eigen::Vector3d({hx_double,hy_double,hz_double}));
}

std::pair<Eigen::Vector3d,Eigen::Vector3d> analytic_value_old(SpaceTimePoint p, double sigma, double eps, double mu, double freq)
{
	auto x = p[0]; auto y = p[1]; auto z = p[2]; auto t = p[3];
	double c = 1/sqrt(eps*mu);
	// double ax=5e-2;
	// double az=10e-2;
	// double ay=2.5e-2; //momentarily useless
	double ax=1;
	double az=1;
	double ay=1; //momentarily useless
	double ksi = sigma/eps/2;
	double alph1 = pow(PI*c/ax,2);
	double alph2 = 0.25*pow(sigma/eps,2);
	double alpha;
	bool flag = true;
	if (alph1>alph2)
	{
		flag = true;
		alpha =  sqrt(alph1-alph2);
	}
	else
	{
		flag = false;
		alpha = sqrt(alph2-alph1);
	}
	
	
	// std::cout << "Debug: " << pow(3.141592*c/ax,2) << " " << pow(sigma/eps,2)/4 << std::endl;
	
	// std::cout << "Parameters:" << std::endl << "C = " << c << "  ksi = " << ksi 
											// << "  alpha = " << alpha << std::endl; 
	
	// int32_t lim1 = floor( (t * c - z) / (2 * az) );
	// int32_t lim2 = floor( (t * c + z) / (2 * az) - 1);
	
	// std::cout << "At time t = " << t << " limits are : " << lim1 << "\t" << lim2 << std::endl;
	
	double a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12;
	int32_t j;
	// i=j=0;
	a1=a2=a3=a4=a5=a6=a7=a8=a9=a10=a11=a12=0;
	
	/******************************HX******************************************************/
	j=0;
	while (true)
	{
		double k = (z + 2*az*j) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			a1 += inverse_laplace_transform_ey_old(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	j=0;
	while (true)
	{
		double k = (2*az*j + 2*az - z) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			 a2 -= inverse_laplace_transform_ey_old(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	
	double hx_double = 0;
	double hy_double = 0;
	double hz_double = 0;
	double ex_double = 0;
	double ey_double = sin(PI*x/ax)*(a1+a2+a3+a4);
	double ez_double = 0;
	
	return std::make_pair<Eigen::Vector3d,Eigen::Vector3d>(Eigen::Vector3d({ex_double,ey_double,ez_double}),Eigen::Vector3d({hx_double,hy_double,hz_double}));
}

std::pair<Eigen::Vector3d,Eigen::Vector3d> analytic_value_cyl(SpaceTimePoint p, double sigma, double eps, double mu, double freq)
{
	auto x = p[0]; auto y = p[1]; auto z = p[2]; auto t = p[3];
	double rho   = std::sqrt(std::pow(x,2)+std::pow(y,2));
	double theta;
	
	if (x==0 && y == 0)
		theta = 0;
	else if (x>=0)
		theta = std::asin(y/rho);
	else
		theta = -std::asin(y/rho)+PI;
	
	double k1 = -207.84727995444730376134397098783; 	//HARD CODED: change if geometry, frequency or material parameters change
	double k2 =   87.296645428299512314918079813159;	//HARD CODED: change if geometry, frequency or material parameters change

	double c = 1/sqrt(eps*mu);
	double ra=  0.05;
	double rb=  0.1;
	double az=  0.6;
	double ksi = sigma/eps/2;
	
	double alph1 = std::pow(c*k2,2);
	double alph2 = 0.25*pow(sigma/eps,2);
	
	if (c < 0.5*c0)
		alph1 = std::pow(c*k1,2);
	
	// if (rho > ra)
	// {
		// std::cout << "rho = " << rho << std::endl;
		// std::cout << "c = " << c << std::endl;
		// std::cout << "k = " << k2 << std::endl;
	// }
	
	double alpha;
	bool flag = true;
	if (alph1>alph2)
	{
		flag = true;
		alpha =  sqrt(alph1-alph2);
	}
	else
	{
		flag = false;
		alpha = sqrt(alph2-alph1);
	}
	
	double a1, a2, a3, a4, a5, a6;
	int32_t j;
	// i=j=0;
	a1=a2=a3=a4=a5=a6=0;
	
	/******************************HX******************************************************/
	j=0;
	while (true)
	{
		double k = (z + 2*az*j) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			a1 += inverse_laplace_transform_hx(t,k, alpha, ksi, c, freq, flag);
			a3 += inverse_laplace_transform_hz(t,k, alpha, ksi, c, freq, flag);
			a5 += inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	j=0;
	while (true)
	{
		double k = (2*az*j + 2*az - z) / c;
		// if ((t-k) > 1e-8*k)
		if (t > k + 1e-18)
		{
			a2 -= inverse_laplace_transform_hx(t,k, alpha, ksi, c, freq, flag);
			a4 += inverse_laplace_transform_hz(t,k, alpha, ksi, c, freq, flag);
			a6 += inverse_laplace_transform_ey(t,k, alpha, ksi, c, freq, flag);
		}
		else
			break;
		j++;
	}
	
	Eigen::Matrix3d Mconvert;
	
	Mconvert << cos(theta), -rho*sin(theta), 0,
			 	sin(theta),  rho*cos(theta), 0,
				0,           0,              1;
	
	Eigen::Vector3d h_cylindric(0,0,0), e_cylindric(0,0,0);

	if (c < 0.5*c0)
	{
		e_cylindric = Eigen::Vector3d({gsl_sf_bessel_J1(k1*rho)*(a1+a2),0,k1*c*gsl_sf_bessel_J0(k1*rho)*(a3+a4)});
		h_cylindric = Eigen::Vector3d({0,sqrt(eps/mu)*gsl_sf_bessel_J1(k1*rho)*(a5+a6),0});
	}
	else
	{
		double const_coeff = (k1/k2)*gsl_sf_bessel_J0(k1*ra)/(gsl_sf_bessel_J0(k2*ra)*gsl_sf_bessel_Y0(k2*rb)-gsl_sf_bessel_J0(k2*rb)*gsl_sf_bessel_Y0(k2*ra));
		// std::cout << "Amp2_1 = " <<  const_coeff*gsl_sf_bessel_Y0(k2*rb) << std::endl;
		// std::cout << "Amp2_2 = " << -const_coeff*gsl_sf_bessel_J0(k2*rb) << std::endl;
		e_cylindric=Eigen::Vector3d({(gsl_sf_bessel_Y0(k2*rb)*gsl_sf_bessel_J1(k2*rho)-gsl_sf_bessel_J0(k2*rb)*gsl_sf_bessel_Y1(k2*rho))*(a1+a2),
		                              0,
									  k2*c*(gsl_sf_bessel_Y0(k2*rb)*gsl_sf_bessel_J0(k2*rho)-gsl_sf_bessel_J0(k2*rb)*gsl_sf_bessel_Y0(k2*rho))*(a3+a4)});
		h_cylindric=Eigen::Vector3d({ 0, 
		                              sqrt(eps/mu)*(gsl_sf_bessel_Y0(k2*rb)*gsl_sf_bessel_J1(k2*rho)-gsl_sf_bessel_J0(k2*rb)*gsl_sf_bessel_Y1(k2*rho))*(a5+a6),
									  0});
		e_cylindric *= const_coeff;
		h_cylindric *= const_coeff;
	}

	
	// Eigen::Vector3d h_cartesian = Mconvert*h_cylindric;
	// Eigen::Vector3d e_cartesian = Mconvert*e_cylindric;
	
	return std::make_pair<Eigen::Vector3d,Eigen::Vector3d>(Mconvert*e_cylindric,Mconvert*h_cylindric);
}
