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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "timecounter.h"
#include "sgnint32_t.hpp"
#include "mapped_file.h"
#include "strtot.hpp"



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
typedef std::string												Meshtype;
typedef std::string												OutputMode;
typedef double													Amplitude;
typedef double 													Frequency;
typedef double 													WaveNumber;
typedef double													Duration;
typedef std::array<double,4>									SpaceTimePoint;
typedef std::array<double,3> 									WaveVector;

//lists of allowed string constants
const std::vector<Primitive>							definables		= {"material","source","mesh","bc","simulation","geometry","output"};
const std::vector<Sourcetype>   						sourcetypes   	= { "e", "b", "j" };
const std::vector<Profile>   							profiles   		= { "wave", "gaussian", "dc" };
const std::vector<Direction>    						directions    	= { "x", "y", "z" };
const std::vector<BaseFunction> 						modes 	= { "sin", "cos" };
const std::vector<BoundaryConditionType>				bctypes			= { "pec", "pmc", "pml" };
const std::vector<Meshtype>								meshtypes		= { "tetrahedral", "cartesian", "none"};
const std::vector<Meshtype>								meshers		    = { "netgen", "gmsh", "none"};
const std::vector<OutputMode>							outputmodes		= { "silo", "probepoint"};

const std::runtime_error pml_missing(std::string("Sorry, PML not implemented yet, getting there!"));
const std::runtime_error pmc_missing(std::string("Sorry, PMC not implemented yet, getting there!"));
const std::runtime_error bc_unknown_type(std::string("Unrecognized boundary condition type! Available: pec, pmc, pml"));
const std::runtime_error bc_unknown_parameter(std::string("Unrecognized boundary condition parameter! Available: type"));	
const std::runtime_error src_unknown_direction(std::string("Unrecognized direction! Available: x, y, z"));
const std::runtime_error src_unknown_type(std::string("Unrecognized source type!") + 
                                          std::string("Available: e (electric field), b (magnetic field), j (current density)"));
const std::runtime_error src_unknown_profile(std::string("Unrecognized source profile! Available: dc, wave, gaussian"));
const std::runtime_error coordinates_syntax(std::string("coordinates must be inside braces {...,..,..}"));
const std::runtime_error unbalanced_bracket(std::string("unbalanced bracket"));
const std::runtime_error too_many_coords(std::string("Maximum of three coordinates!"));
const std::runtime_error too_few_coords(std::string("Undefined end to list of coordinates!"));
const std::runtime_error src_unknown_bf(std::string("Unrecognized base function! Available: sin, cos"));
const std::runtime_error src_unknown_parameter(std::string("unrecognized parameter for source!"));
const std::runtime_error material_unknown_parameter(std::string("Unrecognized material parameter! Available: epsilon, mu, sigma, chi"));
const std::runtime_error mesh_unknown_type(std::string("undefined mesh type! Available: tetrahedral, cartesian"));
const std::runtime_error mesh_unknown_mesher(std::string("undefined mesher! Available: netgen, gmsh"));
const std::runtime_error mesh_unknown_parameter(std::string("undefined mesh parameter! Available: file, name, type"));
const std::runtime_error sim_unknown_output(std::string("undefined output mode type! Available: silo, probe"));
const std::runtime_error sim_unknown_parameter(std::string("undefined simulation parameter! Available: source, mesh, duration, output"));
const std::runtime_error out_unknown_parameter(std::string("undefined output parameter! Available: mode, period, probe, name"));
const std::runtime_error set_wo_define(std::string("define something before setting variables"));
const std::runtime_error unknown_define(std::string("can only define material, mesh, boundary condition or source"));
const std::runtime_error end_wo_define(std::string("ending non defined definition"));
const std::runtime_error unknown_instruction(std::string("Unknown instruction inside define block"));
const std::runtime_error unexpected_end(std::string("File ended unexpectedly before ending definition"));
const std::runtime_error out_of_bounds_freq(std::string("Output frequency must be a value between 0 and 1!"));

template<typename T>
void sort_unique(std::vector<T>& v) //useful as stand-alone
{
	std::sort(v.begin(), v.end());
	auto uniq_iter = std::unique(v.begin(), v.end());
	v.erase(uniq_iter, v.end());
}

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