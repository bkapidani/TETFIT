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
 *       documentation and/or other Materials provided with the distribution.
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

//file Discretization.hpp
#include <iostream>
#include <array>
#include <stdexcept>
#include <stdio.h>
#include <cstdlib>
#include <utility> 
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <map>
#include <iomanip>
#include <set>
// #include <limits>
#include <sstream>
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
const std::vector<Primitive>							definables		= {"material","source","mesh","bc","simulation","geometry"};
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
const std::runtime_error mesh_unknown_parameter(std::string("undefined mesh parameter! Available: file, type"));
const std::runtime_error sim_unknown_output(std::string("undefined output mode type! Available: silo, probe"));
const std::runtime_error sim_unknown_parameter(std::string("undefined simulation parameter! Available: source, mesh, duration"));
const std::runtime_error set_wo_define(std::string("define something before setting variables"));
const std::runtime_error unknown_define(std::string("can only define material, mesh, boundary condition or source"));
const std::runtime_error end_wo_define(std::string("ending non defined definition"));
const std::runtime_error unknown_instruction(std::string("Unknown instruction inside define block"));
const std::runtime_error unexpected_end(std::string("File ended unexpectedly before ending definition"));

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

class BoundaryCondition
{
	public:
	BoundaryCondition()
	: type("none")
	{
		thickness=0; //just for PML, otherwise unused
		is_set = false;
	}
	
	bool Set(void) { return is_set; };
	
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		is_set = true;
		if (param == "type")
		{
			// Temporary error messages
			if (value == "pml")
				MyThrow(input_line,pml_missing);
			else if (value == "pmc")
			{
				type = value;
				this->val=0;
			}
			else if (value == "pec")
			{
				type = value;
				this->val=0;
			}
			else
				MyThrow(input_line,bc_unknown_type);
		}
		else if (param == "thickness")
			thickness = std::stod(value);
		else	
			MyThrow(input_line,bc_unknown_parameter);
	}
	const std::string& 		Type(void) { return type; }
	const double& 			GetThickness(void) { return thickness; }
	const double&			GetValue(void) { return val; }
	
	private:
	std::string type;
	bool is_set;
	double thickness, val;
};

class Source
{
	public:
	Source()
	{
		is_set  = false;
		prof    = "dc";
		dir   	= "x";
		amp   	= 0;
		freq  	= 0;
		width   = 0;
		kvec[0] = kvec[1] = kvec[2] = 0;
		st      = "j";
		bfuncs[0] = bfuncs[1] = bfuncs[2] = "cos";
		center_coords[0] = center_coords[1] = center_coords[2] = 0;
		surface_label   = 0;
	}
	
	
	bool Set(void) { return is_set; };
	
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		is_set = true;
		if (param == "direction")
		{
			if (std::find(directions.begin(),directions.end(),value) == directions.end())
				MyThrow(input_line,src_unknown_direction);
			else
				this->dir = value;
		}
		else if (param == "type")
		{
			if (std::find(sourcetypes.begin(),sourcetypes.end(),value) == sourcetypes.end())
				MyThrow(input_line,src_unknown_type);
			else
				this->st = value;
		}
		else if (param == "profile")
		{
			if (std::find(profiles.begin(),profiles.end(),value) == profiles.end())
				MyThrow(input_line,src_unknown_profile);
			else
				this->prof = value;
		}
		else if (param == "center")
		{
			auto  i = value.begin();
			if (*i != '{')
				MyThrow(input_line,coordinates_syntax);
			else
			{
				uint8_t k=0;
				i++;
				
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					
					while (*i != ',' && *i != '}' && i != value.end())
					{
						coord.push_back(*i);
						i++;
					}
					
					// std::cout << value << std::endl;
					
					if (i == value.end())
						MyThrow(input_line,unbalanced_bracket);
					else 
					{
						if (k < 3)
						{
							center_coords[k]= std::stod(coord);
							k++;
						}
						else
							MyThrow(input_line,too_many_coords);
						
						if (*i == ',')
							i++;
					}
				}
				
				if (i == value.end())
					MyThrow(input_line,too_few_coords);
			}
		}
		else if (param == "wavevector")
		{
			auto  i = value.begin();
			if (*i != '{')
				MyThrow(input_line,coordinates_syntax);
			else
			{
				uint8_t k=0;
				i++;
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					while (*i != ',' && *i != '}' && i != value.end())
					{
						coord.push_back(*i);
						i++;
					}
					if (i == value.end())
						MyThrow(input_line,unbalanced_bracket);
					else 
					{
						if (k < 3)
						{
							kvec[k]= std::stod(coord);
							k++;
						}
						else
							MyThrow(input_line,too_many_coords);
						
						if (*i == ',')
							i++;
					}
					
				}
				
				if (i == value.end())
					MyThrow(input_line,too_few_coords);
			}
		}
		else if (param == "mode")
		{
			auto  i = value.begin();
			if (*i != '{')
				MyThrow(input_line,coordinates_syntax);
			else
			{
				uint8_t k=0;
				i++;
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					while (*i != ',' && *i != '}' && i != value.end())
					{
						if (*i != ' ')
							coord.push_back(*i);
						i++;
					}
					if (i == value.end())
						MyThrow(input_line,unbalanced_bracket);
					else 
					{
						if (k < 3)
						{
							
							if (std::find(modes.begin(),modes.end(),coord) == modes.end())
								MyThrow(input_line,src_unknown_bf);
							
							bfuncs[k] = coord;
							k++;
						}
						else
							MyThrow(input_line,too_many_coords);
						
						if (*i == ',')
							i++;
					}
				}
				
				if (i == value.end())
					MyThrow(input_line,too_few_coords);
			}
		}
		else if (param == "Xmode")
		{
			if (std::find(modes.begin(),modes.end(),value) == modes.end())
				MyThrow(input_line,src_unknown_bf);
			else
				this->bfuncs[0] = value;
		}
		else if (param == "Ymode")
		{
			if (std::find(modes.begin(),modes.end(),value) == modes.end())
				MyThrow(input_line,src_unknown_bf);
			else
				this->bfuncs[1] = value;
		}
		else if (param == "Zmode")
		{
			if (std::find(modes.begin(),modes.end(),value) == modes.end())
				MyThrow(input_line,src_unknown_bf);
			else
				this->bfuncs[2] = value;
		}
		else if (param == "amplitude")
			this->amp = std::stod(value);
		else if (param == "frequency")
			this->freq = std::stod(value);
		else if (param == "width")
			this->width = std::stod(value);
		else if (param == "kx")
			this->kvec[0] = std::stod(value);
		else if (param == "ky")
			this->kvec[1] = std::stod(value);
		else if (param == "kz")
			this->kvec[2] = std::stod(value);
		else if (param == "surface")
			this->surface_label = std::stod(value);
		else
			MyThrow(input_line,src_unknown_parameter);
	}
	
	const Sourcetype& 	Type(void) { return st; }
	const uint32_t&     Surface(void) { return surface_label; } 
	const Direction& 	GetDirection(void) { return dir; }
	const WaveNumber&	Getkx(void) { return kvec[0]; }
	const WaveNumber&	Getky(void) { return kvec[1]; }
	const WaveNumber&	Getkz(void) { return kvec[2]; }
	const Frequency&	GetFreq(void) { return freq; }
	const Amplitude&	GetAmp(void) { return amp; }
	
	Eigen::Vector3d Compute(SpaceTimePoint p) //returns a vector, so one can use superposition of sources
	{
		double ret = amp*cos(2*PI*freq*p[3]); // if the source is dc, we are already done!
		
		if (prof == "gaussian")
		{
			double exponent = - ((pow(p[0]-center_coords[0],2)+
			                      pow(p[1]-center_coords[1],2)+
							      pow(p[2]-center_coords[2],2))/
								  (2*pow(width,2)));
			ret *= exp(exponent);
		}
		else if (prof == "wave")
		{
			// std::cout << "Fai la cosa giusta" << std::endl;
			
			for (uint8_t j=0; j<3; j++)
			{
				if (bfuncs[j] == "sin")
					ret *= sin(2*PI*kvec[j]*(p[j]-center_coords[j]));
				else
					ret *= cos(2*PI*kvec[j]*(p[j]-center_coords[j]));
			}
		}
		
		if (dir == "x")
			return Eigen::Vector3d({ret,0,0});
		else if (dir == "y")
			return Eigen::Vector3d({0,ret,0});
		else
			return Eigen::Vector3d({0,0,ret});
	}
	
	private:
	Sourcetype st;
	Profile prof;
	Frequency freq;
	Amplitude amp;
	double    width;
	Direction dir;
	WaveVector kvec, center_coords;
	std::array<BaseFunction,3> bfuncs;
	uint32_t surface_label;
	bool is_set;
	// double (*foo)();
};

class Material
{
	public:
	Material()
	{
		epsilon = epsilon0;
		sigma = 0;
		mu = mu0;
		chi = 0;
	}
	
	//setters
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "epsilon")
			epsilon = epsilon0*std::stod(value);
		else if (param == "mu")
			mu = mu0*std::stod(value);
		else if (param == "sigma")
			sigma = std::stod(value);
		else if (param == "chi")
			chi = std::stod(value);
		else
			MyThrow(input_line,material_unknown_parameter);
	}

	//getters
	double Epsilon(void) { return epsilon; }
	double Mu(void) { return mu; }
	double Sigma(void) { return sigma; }
	double Chi(void) { return chi; }
	
	private:
	double epsilon, sigma, mu, chi;
};

class Mesh
{
	public:
	Mesh() 
	: type("none"), mesher("none")
	{
		loaded = false;
		xstep = ystep = zstep = 0;
	}

	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "file")
			file = value;
		else if (param == "type")
		{
			if (std::find(meshtypes.begin(),meshtypes.end(),value) == meshtypes.end())
				MyThrow(input_line,mesh_unknown_type);
			type = value;
		}
		else if (param == "mesher")
		{
			if (std::find(meshers.begin(),meshers.end(),value) == meshers.end())
				MyThrow(input_line,mesh_unknown_mesher);
			mesher = value;
		}
		else if (param == "xstep")
			xstep = std::stod(value);
		else if (param == "ystep")
			ystep = std::stod(value);
		else if (param == "zstep")
			zstep = std::stod(value);
		else
			MyThrow(input_line,mesh_unknown_parameter);
	}
	
	const std::string& GetFileName() { return file; }
	const std::string& GetMeshType() { return type; }
	const std::string& GetMesher()   { return mesher; }
	const double& GetLx() { return xstep; }
	const double& GetLy() { return xstep; }
	const double& GetLz() { return xstep; }
	bool IsLoaded() { return loaded; }
	void Switch() { loaded = !loaded; }
	
	private:
	std::string file;
	std::string type;
	std::string mesher;
	bool loaded;
	double xstep,ystep,zstep; //used only when mesh type is cartesian
};

class Simulation
{
	public:
	Simulation()
	{
		d=0;
		sources=std::vector<uint32_t>({1});
		mesh_label=1;
		mode = "silo";
	}
	
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "source")
		{
			sources.push_back(std::stod(value));
			sort_unique(sources);
		}
		else if (param == "mesh")
			mesh_label = std::stod(value);
		else if (param == "duration")
			d = std::stod(value);
		else if (param == "output")
		{
			if (std::find(outputmodes.begin(),outputmodes.end(),value) == outputmodes.end())
				MyThrow(input_line,sim_unknown_output);
			mode = value;
		}
		else if (param == "probe")
		{
			auto  i = value.begin();
			if (*i != '{')
				MyThrow(input_line,coordinates_syntax);
			else
			{
				uint8_t k=0;
				i++;
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					while (*i != ',' && *i != '}' && i != value.end())
					{
						coord.push_back(*i);
						i++;
					}
					if (i == value.end())
						MyThrow(input_line,unbalanced_bracket);
					else 
					{
						if (k < 3)
						{
							probepoint[k]= std::stod(coord);
							k++;
						}
						else
							MyThrow(input_line,too_many_coords);
						
						if (*i == ',')
							i++;
					}
					
				}
				
				if (i == value.end())
					MyThrow(input_line,too_few_coords);
			}
		}
		else
			MyThrow(input_line,sim_unknown_parameter);
	}
	
	const Duration& Time(void) const { return d; } 
	const uint32_t& MeshLabel(void) const { return mesh_label; }
	const OutputMode& Output(void) const { return mode; }
	const Eigen::Vector3d& Probe(void) const { return probepoint; }
	
	private:
	Duration d;
	Eigen::Vector3d probepoint;
	std::vector<uint32_t> sources; //can combine multiple sources
	uint32_t mesh_label;
	OutputMode mode;
};

class Discretization
{	
	public:
	Discretization(std::string inputfile)
	{
		std::string line;
		input_line = 1;
		std::ifstream ReadFile;//(inputfile.c_str());
		ReadFile.open(inputfile.c_str());
		bool in_definition = false;
		bool mesh_added = false;
		
		
		char action[10], token[20], value[64];
		std::string thing_being_defined;
		uint32_t definition_label,input_line;
		
		// std::cout << "ciao!" << std::endl;
		input_line = 1;
		while(getline(ReadFile,line))
		{
			// std::cout << "ciao!" << std::endl;
			auto c_line = line.c_str();
			if (line.size()>0    &&
				c_line[0] != '#' &&       /* lines that begin with '#' are comments         */
				c_line[0] != '\n'  )      /* empty lines for better readability are allowed */
			{
				sscanf(c_line,"%s %s %s",action,token,value);
				std::string instr(action), tok(token), val(value);
				
				if (instr[0] != '#')
				{
					
					if (!in_definition)
					{
						if (instr != "DEFINE")
							MyThrow(input_line,set_wo_define);
						else 
						{
							if (std::find(definables.begin(),definables.end(),tok) == definables.end())
								MyThrow(input_line,unknown_define);
							else
							{
								
								thing_being_defined = tok;
								definition_label = std::stod(val);
								in_definition = true;
								
								if (thing_being_defined == "material")
									Materials[definition_label]=Material();
								else if (thing_being_defined == "bc")
									BCs[definition_label]=BoundaryCondition();
								else if (thing_being_defined == "mesh")
									Meshes[definition_label]=Mesh();
								else if (thing_being_defined == "source")
									Sources[definition_label]=Source();
								else if (thing_being_defined == "simulation")
									Simulations[definition_label]=Simulation();
							}
						}
					}
					else if (instr == "END")
					{
						if (tok != thing_being_defined || std::atoi(value) != definition_label)
							MyThrow(input_line,end_wo_define);
						else 
							in_definition = false;
					}
					else if (instr != "SET")
						MyThrow(input_line,unknown_instruction);
					else
					{
						// std::cout << "We have: " << std::endl;
						// std::cout << '\t' << Meshes.size()      << std::setw(30) << " meshes"              << std::endl;
						// std::cout << '\t' << Materials.size()   << std::setw(30) << " materials"           << std::endl;
						// std::cout << '\t' << BCs.size()         << std::setw(30) << " boundary conditions" << std::endl;
						// std::cout << '\t' << Sources.size()     << std::setw(30) << " sources"             << std::endl;
						// std::cout << '\t' << Simulations.size() << std::setw(30) << " simulations"         << std::endl;
						
						if (value[0] == '{') //check for spaces
						{
							auto i = line.begin();
							val.clear();
							
							while (*i != '{')
								i++;
							while (*i != '}' && i != line.end())
							{
								val.push_back(*i);
								i++;
							}
							
							if (i == line.end())
								MyThrow(input_line,unbalanced_bracket);
							else
								val.push_back(*i);
						}
						
						if (thing_being_defined == "material")
							Materials[definition_label].SetParam(input_line,tok,val);
						else if (thing_being_defined == "bc")
							BCs[definition_label].SetParam(input_line,tok,val);
						else if (thing_being_defined == "mesh")
							Meshes[definition_label].SetParam(input_line,tok,val);
						else if (thing_being_defined == "source")
							Sources[definition_label].SetParam(input_line,tok,val);
						else if (thing_being_defined == "simulation")
							Simulations[definition_label].SetParam(input_line,tok,val);
						else
							MyThrow(input_line,unknown_define);
					}
				
				}
				
				// std::cout << std::setw(20) << instr << std::setw(20) << tok << std::setw(20) << val << std::endl;
			}
			
			input_line++;
		}
		
		if (in_definition)
			MyThrow(input_line,unexpected_end);
	
		ReadFile.close();
	}

	void Run(void)
	{
		timecounter t_sim;
		t_sim.tic();
		
		for (auto sims : Simulations)
			RunSimulation(sims.second);
		t_sim.toc();
		
		std::cout << std::endl;
		std::cout << "----------- All simulations ran. Elapsed time: " << t_sim << " seconds ----------" << std::endl;
	}
	
	void FlushMesh(void)
	{
		auto m = Meshes[loaded_mesh_label];
		if (m.IsLoaded())
		{
			if (m.GetMeshType() == "tetrahedral")
			{
				// std::cout << "Ma qui?" << std::endl;
				std::vector<cluster_list>().swap(vtf_list);
				std::vector<cluster_list>().swap(ftv_list);
				std::vector<cluster_list>().swap(fte_list);
				std::vector<cluster_list>().swap(etf_list);
				std::vector<cluster_list>().swap(etn_list);
				std::vector<cluster_list>().swap(nte_list);
				std::vector<uint32_t>().swap(vol_material);
				std::vector<std::vector<uint32_t>>().swap(associated_volumes);
				std::vector<uint8_t>().swap(classify_edges);
				std::vector<uint8_t>().swap(classify_surfaces);
				std::vector<volume_type>().swap(volumes);
				std::vector<surface_type>().swap(surfaces);
				std::vector<edge_type>().swap(edges);
				std::vector<Eigen::Vector3d>().swap(pts);
			}
			else
			{
				//Still need to merge variables
			}
			m.Switch();
		}
	}
	
	void RunSimulation(const Simulation& s)
	{
		auto m = Meshes[s.MeshLabel()];
		Duration simulation_time;
		
		meshlock.lock(); //lock the access to the meshes map
		if (!m.IsLoaded())
		{
			FlushMesh();
			ReadMesh(m);
			t_step = (double(9)/double(20))*estimate_time_step_bound();
			
			ConstructMaterialMatrices();
		}
		
		// std::ofstream debug_file("massmatrix.dat");
		
		// for (auto s : surfaces)
			// debug_file << "Triangle: [ " << std::get<0>(s) << " " << std::get<1>(s) << " " << std::get<2>(s) << " ]" << std::endl;
		// for (auto e : edges)
			// debug_file << "Edge: [ "<< std::get<0>(e) << " " << std::get<1>(e) << " ]" << std::endl;
		
		// debug_file.close();
		
		for (auto bb : BCs)
			std::cout << "Boundary condition: " << bb.first << " of type " << bb.second.Type() << std::endl;
		for (auto bb : Sources)
			std::cout << "Source: " << bb.second.Surface() << " of type " << bb.second.Type() << std::endl;
		simulation_time = s.Time();
		
		std::cout  << std::endl << std::endl; 
		std::cout <<"------------------------ Running FDTD simulation ------------------------" << std::endl << std::endl;
		double step_time_average=0;
		const uint32_t N_of_steps=simulation_time/t_step;
		size_t i;

		std::cout << std::setw(20) << "Mesh: "            << std::setw(20) << m.GetFileName()              << std::endl;
		std::cout << std::setw(20) << "Simulation time: " << std::setw(20) << simulation_time              << " seconds" << std::endl;
		std::cout << std::setw(20) << "Time step: "       << std::setw(20) << t_step                       << " seconds" << std::endl;
		std::cout << std::setw(20) << "Unknowns: "        << std::setw(20) << U_frac_size+F_frac_size      << std::endl << std::endl;
		
		
		// T time_function;
		timecounter step_cost;
		std::vector<double> numeric_values,numeric_times;
		U = I = Eigen::VectorXd::Zero(edges_size());
		F = Eigen::VectorXd::Zero(surfaces_size());
		Eigen::VectorXd curl_u(surfaces_size()), curl_f(edges_size());
		Eigen::VectorXd U_frac = Eigen::VectorXd::Zero(U_frac_size);
		Eigen::VectorXd F_frac = Eigen::VectorXd::Zero(F_frac_size);
		
		auto start_of_u = U_frac.data();
		auto start_of_f = F_frac.data();
		
		Eigen::Map<Eigen::VectorXd> U_a(start_of_u,H_size), U_b(start_of_u+H_size,P_size), U_c(start_of_u+H_size+P_size,Q_size);
		Eigen::Map<Eigen::VectorXd> F_a(start_of_f,N_size), F_b(start_of_f+N_size,R_size), F_c(start_of_f+N_size+R_size,S_size);
		
		auto full_H = Eigen::MatrixXd(this->H);
		auto full_N = Eigen::MatrixXd(this->N);
		// std::cout << "Max(H) = " << full_H.lpNorm<Eigen::Infinity>() << std::endl;
		// std::cout << "Max(N) = " << full_N.lpNorm<Eigen::Infinity>() << std::endl;
		double current_time;
		
		std::ofstream debug_file("massmatrix.dat");
		
		debug_file << full_H;
		
		debug_file.close();
		
		for (i=0; i*t_step <= simulation_time; i++)
		{
			step_cost.tic();
			
			current_time = double(i)*t_step;
			for (uint32_t ee = 0; ee < edges_size(); ee++)
			{
				if (edge_bcs[ee] != 0  && BCs[edge_bcs[ee]].Type() != "none")
					U[ee] = ComputeEdgeBC(ee,current_time);
				else if (edge_src[ee].size()>0)
				{
					// I[ee] = ComputeCurrentSource(ee,double(i)*t_step);
					U[ee] = ComputeEfieldSource(ee,current_time);
				}
				

			}
			

			// Magnetic Part:
			curl_u     = C*U;
			F_a       -= t_step*N*curl_u;
			F_b        = R*F_frac - t_step*Tr*curl_u;
			F_c        = S*F_frac - t_step*Ts*curl_u;
			F          = T*F_frac;

			// for (uint32_t ff = 0; ff < surfaces_size(); ff++)
			// {
				// if (face_src[ff].size()>0)
					// F[ff] = ComputeBfieldSource(ff,double(i)*t_step);
				// if (face_bcs[ff].size()>0)
					// F[ff] = ComputeFaceBC(ff,double(i)*t_step,face_bcs[ff]);
			// }
			
			//Debug
			std::cout << "Time: "      << std::setw(20) << current_time << '\t'; 
			std::cout << "Maximum F: " << std::setw(20) << F_a.lpNorm<Eigen::Infinity>() << '\t'; 
			std::cout << "Maximum U: " << std::setw(20) << U_a.lpNorm<Eigen::Infinity>() << std::endl;
			
			// Electric Part:			
			curl_f     = C.transpose()*F-I;
			U_a       += t_step*H*curl_f;
			U_b        = P*U_frac + t_step*Mp*curl_f;
			U_c        = Q*U_frac + t_step*Mq*curl_f;
			U          = M*U_frac;


			
			/*auto num_val = GetElectricField(probe_elem);
			numeric_values.push_back(num_val(1));
			numeric_times.push_back(i*t_step);*/
			step_cost.toc();
			step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();

			// ExportFields(s);

			if ((i+1) % 140 == 0)
				std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(9) << step_time_average << "s, " 
						  << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
		}
		
		meshlock.unlock(); //unlock the access to the meshes map

		/* Output stats and fields*/
		// std::ofstream os;
		// os.open("numeric_FIT.dat");
		// for (size_t k=0; k < numeric_values.size(); k++)
		  // os << numeric_times[k] << " " << numeric_values[k] << std::endl;
		// os.close();
		std::cout << std::setw(20) << "Time step average cost is "  << step_time_average/(double(i)) 
		          << " seconds (" << i << " time steps!)" << std::endl;
		std::cout << std::setw(20) << "Total running time is "      << step_time_average << " seconds" << std::endl;
	}
	
	double ComputeEdgeBC(uint32_t e, double t) 
	{ 
		return BCs[edge_bcs[e]].GetValue(); 
	}

	double ComputeEfieldSource(uint32_t e, double t)
	{
		auto eb = edge_barycenter(e);
		SpaceTimePoint p;
		p[0] = eb[0];
		p[1] = eb[1];
		p[2] = eb[2];
		p[3] = t;
		
		Eigen::Vector3d vector_val(0,0,0);
		
		for (auto src : edge_src[e])
		{
			// std::cout << Sources[src].Type() << std::endl;
			if (Sources[src].Type() == "e")
				vector_val += Sources[src].Compute(p);
		}
		
		// std::cout << "-------------" << std::endl << vector_val << std::endl;
		
		return vector_val.dot(pts[abs(etn_list[e][1])]-pts[abs(etn_list[e][0])]);
	}
	// double ComputeFaceSource()
	
	void unique(std::vector<label_surface_type>& arr, std::vector<uint32_t>& new_labels)
	{
		if (!arr.size())
			throw std::invalid_argument("Array must be nonempty!");
		
		uint32_t left = 0;
		uint32_t right = arr.size()-1;

		struct {
			bool operator()(const label_surface_type& t1, const label_surface_type& t2)
			{
				return (t1.first < t2.first);
			}
		} surfcomp;
		
		std::sort(arr.begin(),arr.end(),surfcomp);	
		// std::vector<uint32_t> new_labels(arr.size(),0);

		// std::cout << std::endl;
		
		uint32_t itor = left;
		new_labels[arr[0].second]=itor;
		surfaces.push_back(arr[0].first);
		
		left++;
		while (left<=right)
		{
			if (arr[left].first == arr[left-1].first)
			{
				new_labels[arr[left].second]=itor;
			}
			else
			{
				itor++;
				new_labels[arr[left].second]=itor;
				surfaces.push_back(arr[left].first);
			}
			
			left++;
		}
		
		// arr = new_arr;
		// labels = new_labels;
		
		return;
		
	}
	
	void unique(std::vector<label_edge_type>& arr, std::vector<uint32_t>& new_labels)
	{
		if (!arr.size())
			throw std::invalid_argument("Array must be nonempty!");
		
		uint32_t left = 0;
		uint32_t right = arr.size()-1;

		struct {
			bool operator()(const label_edge_type& t1, const label_edge_type& t2)
			{
				return (t1.first < t2.first);
			}
		} edgecomp;
		
		std::sort(arr.begin(),arr.end(),edgecomp);	
		// std::vector<uint32_t> new_labels(arr.size(),0);
		
		uint32_t itor = left;
		new_labels[arr[0].second]=itor;
		edges.push_back(arr[0].first);
		
		left++;
		while (left<=right)
		{
			if (arr[left].first == arr[left-1].first)
			{
				new_labels[arr[left].second]=itor;
			}
			else
			{
				itor++;
				new_labels[arr[left].second]=itor;
				edges.push_back(arr[left].first);
			}
			
			left++;
		}
		
		// arr = new_arr;
		// labels = new_labels;
		
		return;
		
	}
	
	std::array<Eigen::Vector3d,6> dual_area_vectors(const uint32_t& vol)
	{
		std::array<Eigen::Vector3d,6> ret;
		std::vector<uint32_t> edgs, fcs;
		
		for (auto ff : vtf_list[vol])
		{
			fcs.push_back(abs(ff));
			for (auto ee : fte_list[abs(ff)])
			{
				edgs.push_back(abs(ee));
			}
		}
		sort_unique(edgs);
		std::swap(edgs[1],edgs[2]);
		std::swap(edgs[1],edgs[3]);
		
		// std::cout << "E qui? " << edgs.size() << std::endl;
		std::array<Eigen::Vector3d, 6> ebs;
		for (uint32_t i = 0; i < 6; i++)
			ebs[i] =  edge_barycenter(edgs[i]);


		
		std::array<Eigen::Vector3d, 4> fbs;
		for (uint32_t i = 0; i < 4; i++)
			fbs[i] = face_barycenter(fcs[i]);
		

		
		Eigen::Vector3d vb = vol_barycenter(vol);
		
		/* Area of quadrilateral ABCD is 0.5*|AC x BD|. Orientation of the dual
		 * area vector must be the same of the primal edge vector i.e. it must
		 * satisfy 'dot(PEV, DAV) >= 0' */
		
		double sign = CellVolumes[vol]<0 ? -0.5 : 0.5;
		ret[0] = ((vb-ebs[0]).cross(fbs[2]-fbs[3]))*sign;
		ret[1] = ((vb-ebs[1]).cross(fbs[1]-fbs[3]))*sign;
		ret[2] = ((vb-ebs[2]).cross(fbs[3]-fbs[1]))*sign;
		ret[3] = ((vb-ebs[3]).cross(fbs[1]-fbs[2]))*sign;
		ret[4] = ((vb-ebs[4]).cross(fbs[2]-fbs[0]))*sign;
		ret[5] = ((vb-ebs[5]).cross(fbs[0]-fbs[1]))*sign;
		
		return ret;
	}

	std::array<Eigen::Vector3d, 4> primal_area_vectors(const uint32_t& vol)
	{
		auto nodes = std::vector<uint32_t>({std::get<0>(volumes[vol]),std::get<1>(volumes[vol]),
		                                    std::get<2>(volumes[vol]),std::get<3>(volumes[vol])});
		

		std::array<Eigen::Vector3d, 6> evecs;
		evecs[0]=pts[nodes[1]]-pts[nodes[0]];
		evecs[1]=pts[nodes[2]]-pts[nodes[1]];
		evecs[2]=pts[nodes[2]]-pts[nodes[0]];
		evecs[3]=pts[nodes[3]]-pts[nodes[0]];
		evecs[4]=pts[nodes[3]]-pts[nodes[1]];
		evecs[5]=pts[nodes[3]]-pts[nodes[2]];
		
		std::array<Eigen::Vector3d, 4> ret;
		ret[0] = evecs[1].cross(evecs[4])/2;
		ret[1] = evecs[2].cross(evecs[3])/2;
		ret[2] = evecs[0].cross(evecs[3])/2;
		ret[3] = evecs[0].cross(evecs[2])/2;
		
		return ret;
	}

	bool ReadMesh(Mesh& msh)
	{	
		timecounter tc, tctot;
		
		if (msh.IsLoaded())
			return true;
		else
			msh.Switch();
		
		/* Open file */
		if (msh.GetFileName().size() == 0)
		{
			std::cout << "Invalid mesh file name" << std::endl;
			return false;
		}
		
		uint32_t	lines, linecount;
		
		mapped_file mf(msh.GetFileName());
		
		// std::cout << " * * * Reading NETGEN format mesh * * * " << std::endl;
		
		tctot.tic();
		
		/************************ Read points ************************/
		linecount = 0;
		
		const char *data = mf.mem();
		char *endptr;
		
		lines = strtot<uint32_t>(data, &endptr);
		
		// pts.reserve(lines);
		std::vector<uint32_t> dummy_ass_vols;
		tc.tic();
		while (linecount < lines)
		{
			if ( (linecount%100000) == 0 )
			{
				std::cout << "Reading points: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}

			auto t = parser::read_point_line<double>(endptr, &endptr);
			
			
			Eigen::Vector3d point(std::get<0>(t),std::get<1>(t),std::get<2>(t));
			pts.push_back(point);
			associated_volumes.push_back(dummy_ass_vols);
			
			/* Do something with that point */
			
			linecount++;
		}
		tc.toc();
		
		std::cout << "Reading points: " << linecount;
		std::cout << "/" << lines << " - " << tc << " seconds" << std::endl;
		
		/************************ Read tetrahedra ************************/
		linecount = 0;
		
		lines = strtot<uint32_t>(endptr, &endptr);
		std::vector< tm_tuple > temp_tet;
		temp_tet.reserve(lines);
		
		tc.tic();
		while (linecount < lines)
		{
			if ( (linecount%100000) == 0 )
			{
				std::cout << "Reading tetrahedra: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}
			
			auto t = parser::read_tetrahedron_line<uint32_t>(endptr, &endptr);
			
			//auto t = parser::read_array<uint32_t, 5>(mf.get_line());
			std::vector<uint32_t> vectet(4);
			vectet[0] = std::get<1>(t);
			vectet[1] = std::get<2>(t);
			vectet[2] = std::get<3>(t);
			vectet[3] = std::get<4>(t);
			std::sort(vectet.begin(),vectet.end());
			
			uint32_t       p0(vectet[0]);
			uint32_t       p1(vectet[1]);
			uint32_t       p2(vectet[2]);
			uint32_t       p3(vectet[3]);
			uint32_t       d(std::get<0>(t));
			
			auto tuple = std::make_tuple(volume_type(p0, p1, p2, p3), d);
			temp_tet.push_back( tuple );
			
			linecount++;
		}
		tc.toc();
		
		std::cout << "Reading tetrahedra: " << linecount;
		std::cout << "/" << lines  << " - " << tc << " seconds" << std::endl;
		
		/************************ Sort ************************/
		std::cout << "Sorting data...";
		std::cout.flush();
		
		tc.tic();
		
		/* sort tetrahedra, make unique and move them in geometry */
			
		struct {
			bool operator()(const tm_tuple& t1, const tm_tuple& t2)
			{
				return (std::get<0>(t1) < std::get<0>(t2));
			}
		} mycomp;

		std::sort(temp_tet.begin(), temp_tet.end(), mycomp);
		
		// vtf_list.resize(lines);
		uint32_t tot=0;	

		std::vector<label_surface_type> temp_tri0;
		temp_tri0.resize(4*lines);
		std::vector<int32_t> vol_signs;
		
		for (auto tet : temp_tet)
		{
			auto t = std::get<0>(tet);
			volumes.push_back(t);
			
			uint32_t       p0(std::get<0>(t)); associated_volumes[p0].push_back(tot);
			uint32_t       p1(std::get<1>(t)); associated_volumes[p1].push_back(tot);
			uint32_t       p2(std::get<2>(t)); associated_volumes[p2].push_back(tot);
			uint32_t       p3(std::get<3>(t)); associated_volumes[p3].push_back(tot);

			temp_tri0[tot]           = std::make_pair(surface_type(p0, p1, p2),tot);
			temp_tri0[tot+lines]     = std::make_pair(surface_type(p0, p1, p3),tot+lines);
			temp_tri0[tot+2*lines]   = std::make_pair(surface_type(p0, p2, p3),tot+2*lines);
			temp_tri0[tot+3*lines]   = std::make_pair(surface_type(p1, p2, p3),tot+3*lines);
			tot++;
		
			Eigen::Vector3d v1 = pts[p1] - pts[p0];
			Eigen::Vector3d v2 = pts[p2] - pts[p0];
			Eigen::Vector3d v3 = pts[p3] - pts[p0];		

			// std::cout << v1 << std::endl << std::endl; 
			// std::cout << v2 << std::endl << std::endl;
			// std::cout << v3 << std::endl << std::endl;
			// std::cout << "-----------"   << std::endl;
			
			auto cross_partial = v2.cross(v3);			
			double vol_partial = v1.dot(cross_partial);
			double vol_vol = vol_partial/double(6);
			int32_t sgn  = vol_vol? 1 : -1;
			vol_signs.push_back(sgn);
			CellVolumes.push_back(vol_vol);
			vol_material.push_back(std::get<1>(tet));
		}
		
		//std::cout << volumes_size() << std::endl;
		
		std::vector<tm_tuple>().swap(temp_tet);
		std::vector<uint32_t> labels(4*lines);
		surfaces.reserve(4*lines);
		unique(temp_tri0, labels); //this also fills the surfaces vector
		std::vector<label_surface_type>().swap(temp_tri0);
		
		ftv_list.resize(surfaces.size());
		// vtf_list.reserve(volumes.size());
		
		for (uint32_t k=0; k<lines; k++)
		{
			// std::cout << labels[k] << " " << ftv_list.size() << std::endl;
			sgnint32_t<int32_t> v1(k,-vol_signs[k]);
			sgnint32_t<int32_t> v2(k,vol_signs[k]);
			
			ftv_list[labels[k]].push_back(v1);
			ftv_list[labels[k+lines]].push_back(v2);
			ftv_list[labels[k+2*lines]].push_back(v1);
			ftv_list[labels[k+3*lines]].push_back(v2);
			
			sgnint32_t<int32_t> f1(labels[k],        -vol_signs[k]);
			sgnint32_t<int32_t> f2(labels[k+lines],   vol_signs[k]);
			sgnint32_t<int32_t> f3(labels[k+2*lines],-vol_signs[k]);
			sgnint32_t<int32_t> f4(labels[k+3*lines], vol_signs[k]);
			
			std::vector<sgnint32_t<int32_t>> dummy(4);
			vtf_list.push_back(dummy);
			
			// vtf_list[k][0] = f1; 
			// vtf_list[k][1] = f2; 
			// vtf_list[k][2] = f3;
			// vtf_list[k][3] = f4;
			
			vtf_list[k][0] = f4; 
			vtf_list[k][1] = f3; 
			vtf_list[k][2] = f2;
			vtf_list[k][3] = f1;
		}
		
		std::vector<uint32_t>().swap(labels);
		lines = surfaces.size();
		// intersurface.resize(lines,0);
		// fte_list.resize(lines);
		tot = 0;
		std::vector<label_edge_type> temp_edge0(3*lines);
		
		for (auto t : surfaces)
		{
			uint32_t       p0(std::get<0>(t));
			uint32_t       p1(std::get<1>(t));
			uint32_t       p2(std::get<2>(t));

			temp_edge0[tot]          =  std::make_pair(edge_type(p0, p1),tot);
			temp_edge0[tot+lines]    =  std::make_pair(edge_type(p0, p2),tot+lines);
			temp_edge0[tot+2*lines]  =  std::make_pair(edge_type(p1, p2),tot+2*lines);
			tot++;
		}
		
		std::vector<uint32_t> e_labels(3*lines);
		unique(temp_edge0, e_labels);
		std::vector<label_edge_type>().swap(temp_edge0);
		etf_list.resize(edges.size());
		// fte_list.reserve(surfaces.size());
		// etn_list.resize(edges.size());
		// physical_edges.resize(edges_size(),0);
		
		std::vector<double_triplet> tripletList;
		tripletList.reserve(3*surfaces_size());
		Eigen::SparseMatrix<double> C(surfaces_size(),edges_size());
		n_index.resize(surfaces_size());
		r_index = n_index;
		// s_index.resize(surfaces_size());
		N_size = R_size = S_size = 0;
		std::vector<uint32_t> dual_is_fractured(nodes_size()), primal_is_fractured(volumes_size());
		uint32_t tot_primal_fractured=0;
		uint32_t tot_dual_fractured=0;
		
		for (uint32_t k=0; k<lines; k++)
		{
			sgnint32_t<int32_t> f1(k, 1);
			sgnint32_t<int32_t> f2(k,-1);
			sgnint32_t<int32_t> f3(k, 1);
			
			etf_list[e_labels[k]].push_back(f1);
			etf_list[e_labels[k+lines]].push_back(f2); 
			etf_list[e_labels[k+2*lines]].push_back(f3);

			// std::cout << labels[k] << " " << labels[k+lines] << " " << labels[k+2*lines] << std::endl;
			
			sgnint32_t<int32_t> e1(e_labels[k],1);
			sgnint32_t<int32_t> e2(e_labels[k+lines],-1);
			sgnint32_t<int32_t> e3(e_labels[k+2*lines],1);
			
			std::vector<sgnint32_t<int32_t>> dummy(3);		
			fte_list.push_back(dummy);
			
			fte_list[k][0] = e1;
			fte_list[k][1] = e2;
			fte_list[k][2] = e3;
			
			tripletList.push_back(double_triplet(k,abs(e1), 1));
			tripletList.push_back(double_triplet(k,abs(e2),-1));
			tripletList.push_back(double_triplet(k,abs(e3), 1));
			
			auto vols = ftv_list[k];
			
			
			bool recombine = true;
			switch (vols.size()) 
			{
				case 2: 
				{
					auto vol1= abs(vols[0]);
					auto vol2= abs(vols[1]);
					auto mu_vol1 = Materials[vol_material[vol1]].Mu();
					auto mu_vol2 = Materials[vol_material[vol2]].Mu();
					auto chi_vol1 = Materials[vol_material[vol1]].Chi();
					auto chi_vol2 = Materials[vol_material[vol2]].Chi();
					
					
					if (chi_vol1 != 0 || chi_vol2 != 0)
					{
						//We have magnetic losses
						if (chi_vol1 != chi_vol2)
							recombine=false;
						else if (mu_vol1 != mu_vol2)
							recombine=false;
					}
					
					break;
				}
				case 1:
				{
					break;
				}
				case 0:
				{
					throw std::invalid_argument("Conductor boundary cannot be on mesh boundary!");
					break;
				}
			}
			
			if (recombine && Materials[vol_material[abs(vols[0])]].Chi() == 0)
			{
				classify_surfaces.push_back(1);
				n_index[k]=N_size++;
				
			}
			else if (recombine)
			{
				classify_surfaces.push_back(2);
				r_index[k]=R_size++;
			}
			else
			{
				classify_surfaces.push_back(3);
				
				auto vol1= abs(vols[0]);
				auto n_star1 = abs(vtf_list[vol1]);
				// auto lb1 = std::lower_bound(n_star1.begin(),n_star1.end(),k);
				
				if (!primal_is_fractured[vol1])
				{
					primal_is_fractured[vol1]=S_size+1;
					// s_index[k].push_back(S_size-1+std::distance(n_star1.begin(),lb1));
					S_size+=4;
				}
				// else
					// s_index[k].push_back(primal_is_fractured[vol1]-1 + std::distance(n_star1.begin(),lb1));

				if (vols.size()==2)
				{
					auto vol2= abs(vols[1]);
					auto n_star2 = abs(vtf_list[vol2]);
					
					// auto lb2 = std::lower_bound(n_star2.begin(),n_star2.end(),k);
					
					if (!primal_is_fractured[vol2])
					{
						primal_is_fractured[vol2]=S_size+1;
						// s_index[k].push_back(S_size+std::distance(n_star2.begin(),lb2));
						S_size+=4;
					}
					// else
						// s_index[k].push_back(primal_is_fractured[vol2]-1 + std::distance(n_star2.begin(),lb2));
				
				}
			}
		}
		
		std::vector<uint32_t>().swap(e_labels);
		nte_list.resize(pts.size());
		h_index.resize(edges.size());
		p_index=h_index;
		// q_index.resize(edges.size());
		C.setFromTriplets(tripletList.begin(),tripletList.end());
		this->C=std::move(C);
		
		std::vector<double_triplet>().swap(tripletList);
		
		// std::vector<double_triplet> S_a_trip,S_b_trip,S_c_trip;
		H_size = Q_size = P_size = 0;
		// etn_list.reserve(edges.size());
		
		for (uint32_t k = 0; k < edges.size(); k++)
		{
			auto t = edges[k];
			
			uint32_t       p0(std::get<0>(t));
			uint32_t       p1(std::get<1>(t));
			
			sgnint32_t<int32_t> n1(p0,-1);
			sgnint32_t<int32_t> n2(p1, 1);
			
			std::vector<sgnint32_t<int32_t>> dummy(2);		
			etn_list.push_back(dummy);
			etn_list[k][0] = n1;
			etn_list[k][1] = n2;
			
			sgnint32_t<int32_t> e1(k,-1);
			sgnint32_t<int32_t> e2(k, 1);
			
			nte_list[p0].push_back(e1);
			nte_list[p1].push_back(e2);
			
			bool recombine = true;
			auto itor1 = associated_volumes[p0].begin();
			auto itor2 = associated_volumes[p1].begin();
			auto mat_sigma_label1 = Materials[vol_material[*itor1]].Sigma();
			auto mat_sigma_label2 = Materials[vol_material[*itor2]].Sigma();
			auto mat_eps_label1 = Materials[vol_material[*itor1]].Epsilon();
			auto mat_eps_label2 = Materials[vol_material[*itor2]].Epsilon();
			itor1++;
			
			
			if (mat_sigma_label1 == mat_sigma_label2)
			{
			
				while (itor1 != associated_volumes[p0].end())
				{

					if (Materials[vol_material[*itor1]].Sigma() != mat_sigma_label1)
					{
						recombine=false;
						break;
					}
					else if (mat_sigma_label1 != 0 && Materials[vol_material[*itor1]].Epsilon() != mat_eps_label1)
					{
						recombine=false;
						break;
					}
					
					itor1++;
					// itor2++;
				}
				
				if (recombine)
				{
					while (itor2 != associated_volumes[p1].end())
					{
						/*std::cout << *itor2 <<  std::flush << " " 
						  << vol_material[*itor2] << std::flush << " " 
						  << Materials[vol_material[*itor2]].Sigma() << std::endl;*/
						
						if (Materials[vol_material[*itor2]].Sigma() != mat_sigma_label2)
						{
							recombine=false;
							break;
						}
						else if (mat_sigma_label2 != 0 && Materials[vol_material[*itor2]].Epsilon() != mat_eps_label2)
						{
							recombine=false;
							break;
						}
						
						itor2++;
						// itor2++;
					}
				}				
			}
			else
				recombine=false;
			
			if (recombine && mat_sigma_label1 == 0)
			{
				classify_edges.push_back(1);
				h_index[k]=H_size++;
				// H_size++;
			}
			else if (recombine)
			{
				classify_edges.push_back(2);
				p_index[k]=P_size++;
				// P_size++;
			}
			else
			{
				classify_edges.push_back(3);
				
				auto n_star1 = abs(nte_list[p0]);
				auto n_star2 = abs(nte_list[p1]);
				// auto lb1 = std::lower_bound(n_star1.begin(),n_star1.end(),k);
				// auto lb2 = std::lower_bound(n_star2.begin(),n_star2.end(),k);
				
				if (!dual_is_fractured[p0])
				{
					dual_is_fractured[p0]=Q_size+1;
					// q_index[k][0]=Q_size+std::distance(n_star1.begin(),lb1);
					Q_size+=n_star1.size();
					// nte_list[p0]
				}
				// else
					// q_index[k][0]=dual_is_fractured[p0]-1 + std::distance(n_star1.begin(),lb1);

				if (!dual_is_fractured[p1])
				{
					dual_is_fractured[p1]=Q_size+1;
					// q_index[k][1]=Q_size+std::distance(n_star2.begin(),lb2);
					Q_size+=n_star2.size();
					// nte_list[p0]
				}
				// else
					// q_index[k][1]=dual_is_fractured[p1]-1 + std::distance(n_star2.begin(),lb2);				
			}
		}
		
		this->dual_is_fractured   = std::move(dual_is_fractured);
		this->primal_is_fractured = std::move(primal_is_fractured);
		std::vector<bool> all_false(edges_size(),false);
		this->is_dirichlet = all_false;
		tc.toc();
		
		std::cout << "done - " << tc << " seconds" << std::endl;
		
		/************************ Read boundary surfaces ************************/
		linecount = 0;
		auto num_of_tets=lines;
		lines = strtot<uint32_t>(endptr, &endptr);
		
		face_bcs.resize(surfaces_size());
		edge_bcs.resize(edges_size());
		
		face_src.resize(surfaces_size());
		edge_src.resize(edges_size());
		edge_bids.resize(edges_size(),0);
		uint32_t trinum=0;
		tc.tic();
		while (linecount < lines)
		{
			/*if (ifs.fail())
			{
				std::cout << "Error while reading boundary surfaces" << std::endl;
				return false;
			}*/
			
			if ( (linecount%50000) == 0 )
			{
				std::cout << "Reading triangles: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}
			
			auto t = parser::read_triangle_line<uint32_t>(endptr, &endptr);

			uint32_t       bid( std::get<0>(t) );
			
			std::vector<uint32_t> vecsurf(3);
			vecsurf[0] = std::get<1>(t);
			vecsurf[1] = std::get<2>(t);
			vecsurf[2] = std::get<3>(t);
			std::sort(vecsurf.begin(),vecsurf.end());
			
			uint32_t       p0(vecsurf[0]);
			uint32_t       p1(vecsurf[1]);
			uint32_t       p2(vecsurf[2]);
			surface_type   tri( p0, p1, p2 );
			
			auto itor = std::lower_bound(surfaces.begin(),surfaces.end(),tri);
			uint32_t face_label = std::distance(surfaces.begin(),itor);
			
			// auto src_id = &Sources[bid];
			// auto bc_id  = BCs[bid];
			
			if (BCs[bid].Type() == "pec") // boundary conditions override sources!
			{	
				for (auto ee : fte_list[face_label])
				{
					edge_bcs[abs(ee)]=bid;
					edge_src[abs(ee)].clear();
					is_dirichlet[abs(ee)] = true;
				}
			}
			else
			{				
				for	(auto src_pair : Sources)
				{
					auto src = src_pair.second;
					if (src.Surface() == bid)
					{
						if (src.Type() == "e" || src.Type() == "j")
						{
							for (auto ee : fte_list[face_label])
								if (!is_dirichlet[abs(ee)])
									edge_src[abs(ee)].push_back(src_pair.first);
						}
					}
				}
			}

			
			linecount++;
		}
		
		// for (uint32_t kk=0; kk < edges_size(); kk++)
		// {
			// if (edge_bids[kk] != 0)
			// {
				// for (auto src_pair : Sources)
				// {
					// auto src = src_pair.second;
					// std::cout << kk << "   " << edge_bids[kk] << "   " << src.Surface() << "   " << src.Type() << std::endl;
					
					// if (src.Surface() == edge_bids[kk])
						// edge_src[kk].push_back(src_pair.first);
				// }
			// }
		// }
		
		
		tc.toc();
		
		
		std::cout << "Reading triangles: " << linecount;
		std::cout << "/" << lines  << " - " << tc << " seconds"  << std::endl;

		for (auto src : edge_src)
		{
			sort_unique(src);
			assert(src.size()==1 || src.size() == 0); //ricordati di toglierlo dopo il debugging
		}
		
		
		// for (auto bcs : edge_bcs)
			// sort_unique(bcs);
		
		tctot.toc();
		
		//std::cout << volumes_size() << std::endl;
		
		// std::cout << cyan << "Total time spent in reading mesh: ";
		// std::cout << tctot << " seconds" << nocolor << std::endl;
		
		// for (const auto& nn : vtf_list)
			// assert(nn.size() == 4);
		
		// for (const auto& nn : ftv_list)
			// assert(nn.size() == 2 || nn.size() == 1);
		
		// for (const auto& nn : fte_list)
			// assert(nn.size() == 3);
		
		// for (const auto& nn : etf_list)
		// {
			// assert(nn.size() >= 2);
		// }
		
		// for (const auto& nn : nte_list)
			// assert(nn.size() >= 3);
		// for (const auto& nn : etn_list)
			// assert(nn.size() == 2);	
		
		return true;
	}
	
	double estimate_time_step_bound()
	{
		double ret=1;
		double max_h=1e6;
		
		//std::cout << volumes_size() << std::endl;
		
        for (auto itor = 0; itor < volumes_size(); itor++)
		{
			//std::cout << "E qui?" << std::endl;
			auto vol = fabs(CellVolumes[itor]);
			auto vol_domain = vol_material[itor];
			double eps_vol = Materials[vol_domain].Epsilon();
			double mu_vol = Materials[vol_domain].Mu();
			double c = 1/sqrt(eps_vol*mu_vol);
			auto fareas = primal_area_vectors(itor);
			
			// std::cout << vol << '\t' << vol_domain << '\t' << c << '\t' << std::endl;
			
			for (auto ff : fareas)
			{
				double diameter = 12*vol/ff.norm();
				double dt = diameter/c/2;
				
				if (dt<ret)
					ret=dt;
				if (diameter < max_h)
					max_h = diameter;
			}
		}
		
		// std::cout << "Minimum diameter: " << max_h << std::endl;
		
		return ret;
	}
	
	Eigen::VectorXd power_method_iteration(Eigen::VectorXd b )
	{
		b = b*(1/b.norm());
		auto cb = C*b;
		auto nicb = N*cb;
		auto ctnib = C.transpose()*nicb;
		return H*ctnib;
	}
	
	double estimate_time_step_bound_algebraic(Eigen::VectorXd b )
	{

		b = power_method_iteration(b);
		double lambda = b.norm();
		double tol;
		double lambda_old = lambda;
		
		uint32_t it = 0;
		
		do{
			b = power_method_iteration(b);
			lambda = b.norm();

			// bstar = arma::trans(b_old);
			
			// lambda = dot(b_old,b) / (dot(b_old,b_old));
			
			// std::cout << "Current eigenval: " << lambda << std::endl;
			
			/*parameters for next iteration*/
			tol = fabs(lambda_old-lambda)/lambda_old;
			lambda_old = lambda;
			it++;
		// }while (it<100);
		}while (tol>=1e-3);
		
		return 2/sqrt(lambda);
	}
	
	void ConstructMaterialMatrices(void)
	{
		std::cout << "Assembling constitutive matrices...";
		std::cout.flush();
		timecounter t_material;
		t_material.tic();
		
		std::vector<bool> mu_computed(volumes_size(),false);
		Eigen::MatrixXd local_E, local_S;
		Eigen::Matrix4d local_M, local_Z;
		std::vector<double_triplet> H_trip, Mp_trip, Mq_trip, M_trip, P_trip, Q_trip;
		std::vector<double_triplet> N_trip, Tr_trip, Ts_trip, T_trip, R_trip, S_trip;  
		uint32_t jj,kk;
		
		// auto local_mag_Id = = Eigen::MatrixXd::Identity(4,4);
		
		for (uint32_t nid=0; nid< pts.size(); nid++ )
		{
			// timecounter t_find;
			// t_find.tic();
			auto n_star = nte_list[nid];
			auto local_E_size = n_star.size();
			auto local_Id = Eigen::MatrixXd::Identity(local_E_size,local_E_size);
			std::vector<uint32_t> global_i;
			
			for (auto ee : n_star)
				global_i.push_back(abs(ee));
			
			local_S = local_E = Eigen::MatrixXd::Zero(local_E_size,local_E_size);
			
			auto vols = associated_volumes[nid];
			bool sigma_node = false;
			
			for (auto vv : vols)
			{
				auto vol_nodes = std::vector<uint32_t>({std::get<0>(volumes[vv]),std::get<1>(volumes[vv]),
				                                        std::get<2>(volumes[vv]),std::get<3>(volumes[vv])});
					
				std::vector<uint32_t> edgs_l2g;
				for (auto ff : vtf_list[vv])
					for (auto ee : fte_list[abs(ff)])
						edgs_l2g.push_back(abs(ee));	
				sort_unique(edgs_l2g);
				std::swap(edgs_l2g[1],edgs_l2g[2]);
				std::swap(edgs_l2g[1],edgs_l2g[3]);
				
				// auto vol_edges = edges(vv);
				auto vol_domain = vol_material[vv];
				auto face_vecs = primal_area_vectors(vv);
			
				auto elem_volume = fabs(CellVolumes[vv]);
				
				double eps_vol   = Materials[vol_domain].Epsilon();
				// std::cout << eps_vol << std::endl;
				double sigma_vol = Materials[vol_domain].Sigma();
				double chi_vol	 = Materials[vol_domain].Chi();
				
				if (sigma_vol != 0)
					sigma_node = true;
				// eps_vol = 1;
				
				auto D = vtf_list[vv];
				
				uint32_t e1,e2,e3;
				
				uint32_t k=0;
				for (auto vn : vol_nodes)
				{
					if (vn == nid)
						break;
					else
						k++;
				}
				
				std::vector<uint32_t>::iterator it_gi;
				
				if (k==0)
				{	
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[0]);							
					e1    = std::distance(global_i.begin(),it_gi);
					
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[2]);
					e2    = std::distance(global_i.begin(),it_gi);

					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[3]);						
					e3    = std::distance(global_i.begin(),it_gi);
					
					local_E(e1,e1)=local_E(e1,e1)+(face_vecs[1]*(-2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[1]*
												  (-2.0/3.0/elem_volume*(D[1].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e2)=local_E(e2,e2)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[2]*
												  (-2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e3)=local_E(e3,e3)+(face_vecs[3]*(-2.0/3.0/elem_volume*(D[3].Sgn()))).dot(face_vecs[3]*
												  (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e1)=local_E(e1,e2)=local_E(e1,e2)+(face_vecs[1]*(-2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[2]*
												  (-2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e1)=local_E(e1,e3)=local_E(e1,e3)+(face_vecs[1]*(-2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[3]*
												  (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e2)=local_E(e2,e3)=local_E(e2,e3)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[3]*
												  (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
										
					if (sigma_vol != 0)
					{
						local_S(e1,e1)=local_S(e1,e1)+(face_vecs[1]*(-2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[1]*
													  (-2.0/3.0/elem_volume*(D[1].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e2)=local_S(e2,e2)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[2]*
													  (-2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e3)=local_S(e3,e3)+(face_vecs[3]*(-2.0/3.0/elem_volume*(D[3].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e1)=local_S(e1,e2)=local_S(e1,e2)+(face_vecs[1]*(-2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[2]*
													  (-2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e1)=local_S(e1,e3)=local_S(e1,e3)+(face_vecs[1]*(-2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e2)=local_S(e2,e3)=local_S(e2,e3)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
													  
					}
				}
				else if (k==1)
				{
					
					
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[0]);						
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[1]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[4]);
					e3 = std::distance(global_i.begin(),it_gi);				
					
					
					local_E(e1,e1)=local_E(e1,e1)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[0]*
												  (2.0/3.0/elem_volume*(D[0].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e2)=local_E(e2,e2)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[2]*
												  (-2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e3)=local_E(e3,e3)+(face_vecs[3]*(-2.0/3.0/elem_volume*(D[3].Sgn()))).dot(face_vecs[3]*
												  (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e1)=local_E(e1,e2)=local_E(e1,e2)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[2]*
												  (-2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e1)=local_E(e1,e3)=local_E(e1,e3)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[3]*
												  (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e2)=local_E(e2,e3)=local_E(e2,e3)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[3]*
												  (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
							
					if (sigma_vol != 0)
					{
						local_S(e1,e1)=local_S(e1,e1)+(face_vecs[0]*( 2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[0]*
													  (2.0/3.0/elem_volume*(D[0].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e2)=local_S(e2,e2)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[2]*
													  (-2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e3)=local_S(e3,e3)+(face_vecs[3]*(-2.0/3.0/elem_volume*(D[3].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e1)=local_S(e1,e2)=local_S(e1,e2)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[2]*
													  (-2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e1)=local_S(e1,e3)=local_S(e1,e3)+(face_vecs[0]*( 2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e2)=local_S(e2,e3)=local_S(e2,e3)+(face_vecs[2]*(-2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
					}
				}
				else if (k==2)
				{
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[1]);							
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[2]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[5]);							
					e3 = std::distance(global_i.begin(),it_gi);
					
					local_E(e1,e1)=local_E(e1,e1)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[1]*
												 (2.0/3.0/elem_volume*(D[1].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e2)=local_E(e2,e2)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[0]*
												 (2.0/3.0/elem_volume*(D[0].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e3)=local_E(e3,e3)+(face_vecs[3]*(-2.0/3.0/elem_volume*(D[3].Sgn()))).dot(face_vecs[3]*
												 (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e1)=local_E(e1,e2)=local_E(e1,e2)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[0]*
																 (2.0/3.0/elem_volume*(D[0].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e1)=local_E(e1,e3)=local_E(e1,e3)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[3]*
																 (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e2)=local_E(e2,e3)=local_E(e2,e3)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[3]*
																 (-2.0/3.0/elem_volume*(D[3].Sgn()))*eps_vol)*elem_volume/4;
																 
					if (sigma_vol != 0)
					{
						local_S(e1,e1)=local_S(e1,e1)+(face_vecs[1]*( 2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[1]*
													  (2.0/3.0/elem_volume*(D[1].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e2)=local_S(e2,e2)+(face_vecs[0]*( 2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[0]*
													  (2.0/3.0/elem_volume*(D[0].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e3)=local_S(e3,e3)+(face_vecs[3]*(-2.0/3.0/elem_volume*(D[3].Sgn()))).dot(face_vecs[3]*
													  (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e1)=local_S(e1,e2)=local_S(e1,e2)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[0]*
																	 ( 2.0/3.0/elem_volume*(D[0].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e1)=local_S(e1,e3)=local_S(e1,e3)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[3]*
																	 (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e2)=local_S(e2,e3)=local_S(e2,e3)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[3]*
																	 (-2.0/3.0/elem_volume*(D[3].Sgn()))*sigma_vol)*elem_volume/4;
					}
				}
				else if (k==3)
				{
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[3]);						
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[4]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = std::lower_bound(global_i.begin(),global_i.end(),edgs_l2g[5]);							
					e3 = std::distance(global_i.begin(),it_gi);				
					
					
					local_E(e1,e1)=local_E(e1,e1)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[0]*
												  (2.0/3.0/elem_volume*(D[0].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e2)=local_E(e2,e2)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[1]*
												  (2.0/3.0/elem_volume*(D[1].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e3)=local_E(e3,e3)+(face_vecs[2]*(2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[2]*
												  (2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e2,e1)=local_E(e1,e2)=local_E(e1,e2)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[1]*
												  (2.0/3.0/elem_volume*(D[1].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e1)=local_E(e1,e3)=local_E(e1,e3)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[2]*
												  (2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					local_E(e3,e2)=local_E(e2,e3)=local_E(e2,e3)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[2]*
												  (2.0/3.0/elem_volume*(D[2].Sgn()))*eps_vol)*elem_volume/4;
					
					if (sigma_vol != 0)
					{
						local_S(e1,e1)=local_S(e1,e1)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[0]*
													  (2.0/3.0/elem_volume*(D[0].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e2)=local_S(e2,e2)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[1]*
													  (2.0/3.0/elem_volume*(D[1].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e3)=local_S(e3,e3)+(face_vecs[2]*(2.0/3.0/elem_volume*(D[2].Sgn()))).dot(face_vecs[2]*
													  (2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e2,e1)=local_S(e1,e2)=local_S(e1,e2)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[1]*
													  (2.0/3.0/elem_volume*(D[1].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e1)=local_S(e1,e3)=local_S(e1,e3)+(face_vecs[0]*(2.0/3.0/elem_volume*(D[0].Sgn()))).dot(face_vecs[2]*
													  (2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
						local_S(e3,e2)=local_S(e2,e3)=local_S(e2,e3)+(face_vecs[1]*(2.0/3.0/elem_volume*(D[1].Sgn()))).dot(face_vecs[2]*
													  (2.0/3.0/elem_volume*(D[2].Sgn()))*sigma_vol)*elem_volume/4;
					}
				}
				
				if (!mu_computed[vv])
				{
					
					double mu_vol = Materials[vol_material[vv]].Mu();
					
					// std::cout << mu_vol << std::endl;
					
					auto face_vecs = dual_area_vectors(vv);
					auto fids = vtf_list[vv];
					
					std::vector<uint32_t> abs_fids;
					for (auto ff : fids)
						abs_fids.push_back(abs(ff));
					double coeff=36.0/(fabs(CellVolumes[vv])*(1/mu_vol));
					local_M = local_Z = Eigen::Matrix4d::Zero();
					mu_computed[vv]=true;

					
					local_M(0,0)=((face_vecs[0]).dot(face_vecs[0])+(face_vecs[2]).dot(face_vecs[2])+(face_vecs[3]).dot(face_vecs[3]))*coeff;
					local_M(1,1)=((face_vecs[0]).dot(face_vecs[0])+(face_vecs[1]).dot(face_vecs[1])+(face_vecs[4]).dot(face_vecs[4]))*coeff;
					local_M(2,2)=((face_vecs[1]).dot(face_vecs[1])+(face_vecs[2]).dot(face_vecs[2])+(face_vecs[5]).dot(face_vecs[5]))*coeff;
					local_M(3,3)=((face_vecs[3]).dot(face_vecs[3])+(face_vecs[4]).dot(face_vecs[4])+(face_vecs[5]).dot(face_vecs[5]))*coeff;

					local_M(0,1)=(-(face_vecs[1]).dot(face_vecs[2])-(face_vecs[3]).dot(face_vecs[4]))*coeff;
					local_M(0,2)=(-(face_vecs[0]).dot(face_vecs[1])+(face_vecs[3]).dot(face_vecs[5]))*coeff;
					local_M(0,3)=( (face_vecs[0]).dot(face_vecs[4])+(face_vecs[2]).dot(face_vecs[5]))*coeff;

					local_M(1,0)=local_M(0,1);
					local_M(1,2)=(-(face_vecs[0]).dot(face_vecs[2])-(face_vecs[4]).dot(face_vecs[5]))*coeff;
					local_M(1,3)=( (face_vecs[0]).dot(face_vecs[3])-(face_vecs[1]).dot(face_vecs[5]))*coeff;
					
					local_M(2,0)=local_M(0,2);
					local_M(2,1)=local_M(1,2);
					local_M(2,3)=(-(face_vecs[2]).dot(face_vecs[3])-(face_vecs[1]).dot(face_vecs[4]))*coeff;
					
					local_M(3,0)=local_M(0,3);
					local_M(3,1)=local_M(1,3);
					local_M(3,2)=local_M(2,3);
					
					if (chi_vol != 0)
					{
						double coeff2 = 36.0*chi_vol/(fabs(CellVolumes[vv]));
						
						local_Z(0,0)=((face_vecs[0]).dot(face_vecs[0])+(face_vecs[2]).dot(face_vecs[2])+(face_vecs[3]).dot(face_vecs[3]))*coeff2;
						local_Z(1,1)=((face_vecs[0]).dot(face_vecs[0])+(face_vecs[1]).dot(face_vecs[1])+(face_vecs[4]).dot(face_vecs[4]))*coeff2;
						local_Z(2,2)=((face_vecs[1]).dot(face_vecs[1])+(face_vecs[2]).dot(face_vecs[2])+(face_vecs[5]).dot(face_vecs[5]))*coeff2;
						local_Z(3,3)=((face_vecs[3]).dot(face_vecs[3])+(face_vecs[4]).dot(face_vecs[4])+(face_vecs[5]).dot(face_vecs[5]))*coeff2;

						local_Z(0,1)=(-(face_vecs[1]).dot(face_vecs[2])-(face_vecs[3]).dot(face_vecs[4]))*coeff2;
						local_Z(0,2)=(-(face_vecs[0]).dot(face_vecs[1])+(face_vecs[3]).dot(face_vecs[5]))*coeff2;
						local_Z(0,3)=( (face_vecs[0]).dot(face_vecs[4])+(face_vecs[2]).dot(face_vecs[5]))*coeff2;

						local_Z(1,0)=local_Z(0,1);
						local_Z(1,2)=(-(face_vecs[0]).dot(face_vecs[2])-(face_vecs[4]).dot(face_vecs[5]))*coeff2;
						local_Z(1,3)=( (face_vecs[0]).dot(face_vecs[3])-(face_vecs[1]).dot(face_vecs[5]))*coeff2;
						
						local_Z(2,0)=local_Z(0,2);
						local_Z(2,1)=local_Z(1,2);
						local_Z(2,3)=(-(face_vecs[2]).dot(face_vecs[3])-(face_vecs[1]).dot(face_vecs[4]))*coeff2;
						
						local_Z(3,0)=local_Z(0,3);
						local_Z(3,1)=local_Z(1,3);
						local_Z(3,2)=local_Z(2,3);
						
						local_Z *= 0.5*t_step;
					}
						
					auto local_N = (local_M + local_Z).inverse();
					// auto local_N = (local_M + local_Z).llt().solve(local_mag_Id); 
					auto local_R = local_N*(local_M - local_Z);

					uint32_t offset;
					bool is_frac=false;
					if (primal_is_fractured[nid]>0)
					{
						is_frac=true;
						offset = primal_is_fractured[nid]-1;
					}
					
					jj=kk=0;
						
					for (auto j = abs_fids.begin(); j != abs_fids.end(); j++)
					{
						switch (classify_surfaces[*j])
						{
							case 1 :
							{
								for (auto k = abs_fids.begin(); k != abs_fids.end(); k++)
								{
									if (local_N(jj,kk)!=0)
										N_trip.push_back(double_triplet(n_index[*j],*k,local_N(jj,kk)));
									kk++;
								}
								
								// std::cout << "(" << *j << "," << n_index[*j] << ") ";
								std::cout.flush();
								T_trip.push_back(double_triplet(*j,n_index[*j],double(1)));
								break;
							}
							case 2 :
							{
								for (auto k = abs_fids.begin(); k != abs_fids.end(); k++)
								{
									if (local_R(jj,kk)!=0)
									{
										if (classify_surfaces[*k]==1)
											R_trip.push_back(double_triplet(r_index[*j],n_index[*k],local_R(jj,kk)));
										else if (classify_surfaces[*k]==2)
											R_trip.push_back(double_triplet(r_index[*j],N_size+r_index[*k],local_R(jj,kk)));
										else
											R_trip.push_back(double_triplet(r_index[*j],N_size+R_size+offset+kk,local_R(jj,kk)));
									}

									if (local_N(jj,kk)!=0)
										Tr_trip.push_back(double_triplet(r_index[*j],*k,local_N(jj,kk)));
									kk++;
								}
								
								T_trip.push_back(double_triplet(*j,N_size+r_index[*j],double(1)));
								break;
							}
							default :
							{
								T_trip.push_back(double_triplet(*j,N_size+R_size+offset+jj,double(1)));
								break;
							}
						}
						
						if (is_frac)
						{
							for (kk = 0; kk < 4; kk++)
							{
								if (local_R(jj,kk) != 0)
									S_trip.push_back(double_triplet(offset+jj,offset+kk,local_R(jj,kk)));
							
								if (local_N(jj,kk) != 0)
									Ts_trip.push_back(double_triplet(offset+jj,global_i[kk],local_N(jj,kk)));
							}
						}
						
						jj++;
						kk=0;
					}
				}
				
			}
			
			if (sigma_node != 0)
				local_S *= 0.5*t_step;

			jj=kk=0;
			// if (nid==0)
				// std::cout << std::endl << local_E << std::endl << std::endl;
			/* Handling pec edges */
			for (auto j = global_i.begin(); j != global_i.end(); j++)
			{	
				if (is_dirichlet[*j])
				{
					// std::cout << "ci entri o no?" << std::endl;
					local_E.row(jj).setZero();
					local_E.col(jj).setZero();
					local_E.coeffRef(jj,jj)=1; //
					
					// if (nid==0)
						// std::cout << std::endl << local_E << std::endl << std::endl;
					
					if (sigma_node != 0)
					{
						local_S.row(jj).setZero();
						local_S.col(jj).setZero();
					}
				}
				jj++;
			}
			
			// std::cout << std::endl << local_Id << std::endl << std::endl;
			
			auto local_H = (local_E+local_S).inverse(); //faster than inverse
			auto local_P = local_H*(local_E - local_S);
	
			std::cout << std::endl << local_E << std::endl << std::endl;
			// std::cout << std::endl << local_H << std::endl << std::endl;
			
			uint32_t offset;
			bool is_frac=false;
			if (dual_is_fractured[nid]>0)
			{
				is_frac=true;
				offset = dual_is_fractured[nid]-1;
			}
			jj=kk=0;
			
			for (auto j = global_i.begin(); j != global_i.end(); j++)
			{
				if (!is_dirichlet[*j])
				{
					switch (classify_edges[*j])
					{
						case 1 :
						{
							for (auto k = global_i.begin(); k != global_i.end(); k++)
							{
								
								if (!is_dirichlet[*k] && local_H(jj,kk)!=0)
									H_trip.push_back(double_triplet(h_index[*j],*k,local_H(jj,kk)));
								kk++;
							}
							
							M_trip.push_back(double_triplet(*j,h_index[*j],double(1)));
							break;
						}
						case 2 :
						{
							for (auto k = global_i.begin(); k != global_i.end(); k++)
							{
								if (!is_dirichlet[*k] && local_P(jj,kk)!=0)
								{
									if (classify_edges[*k]==1)
										P_trip.push_back(double_triplet(p_index[*j],h_index[*k],local_P(jj,kk)));
									else if (classify_edges[*k]==2)
										P_trip.push_back(double_triplet(p_index[*j],H_size+p_index[*k],local_P(jj,kk)));
									else 
										P_trip.push_back(double_triplet(p_index[*j],H_size+P_size+offset+kk,local_P(jj,kk)));
								}

								if (local_H(jj,kk)!=0)
									Mp_trip.push_back(double_triplet(p_index[*j],*k,local_H(jj,kk)));
								
								kk++;
							}
							
							M_trip.push_back(double_triplet(*j,H_size+p_index[*j],double(1)));
							break;
						}
						default :
						{
							M_trip.push_back(double_triplet(*j,H_size+P_size+offset+jj,double(1)));
							
							break;
						}
					}
				}
				
				if (is_frac)
				{
					for (kk = 0; kk < local_E_size; kk++)
					{
						if (local_P(jj,kk) != 0)
							Q_trip.push_back(double_triplet(offset+jj,offset+kk,local_P(jj,kk)));
					
						if (local_H(jj,kk) != 0)
							Mq_trip.push_back(double_triplet(offset+jj,global_i[kk],local_H(jj,kk)));
					}
					
					
				}
					
				jj++;
				kk=0;
			}
			
			// t_find.toc();
		}
		
		U_frac_size = H_size + P_size + Q_size;
		F_frac_size = N_size + R_size + S_size;
		
		Eigen::SparseMatrix<double> H(H_size,edges_size()), N(N_size,surfaces_size());
		Eigen::SparseMatrix<double> M(edges_size(),U_frac_size), Mq(Q_size,edges_size()), Mp(P_size,edges_size());
		Eigen::SparseMatrix<double> T(surfaces_size(),F_frac_size),Tr(R_size,surfaces_size()),Ts(S_size,surfaces_size());
		Eigen::SparseMatrix<double> P(P_size,U_frac_size),Q(Q_size,U_frac_size);
		Eigen::SparseMatrix<double> R(R_size,F_frac_size),S(S_size,F_frac_size);
			
		add_to_sparse ass;
		overwrite_to_sparse oss;
		
		T.setFromTriplets(T_trip.begin(),T_trip.end(), oss);
		M.setFromTriplets(M_trip.begin(),M_trip.end(), oss);
		
		N.setFromTriplets(N_trip.begin(),N_trip.end(), ass);
		H.setFromTriplets(H_trip.begin(),H_trip.end(), ass);
		Mp.setFromTriplets(Mp_trip.begin(),Mp_trip.end(), ass);
		Mq.setFromTriplets(Mq_trip.begin(),Mq_trip.end(), ass);
		N.setFromTriplets(N_trip.begin(),N_trip.end(), ass);
		P.setFromTriplets(P_trip.begin(),P_trip.end(), ass);
		Q.setFromTriplets(Q_trip.begin(),Q_trip.end(), ass);
		R.setFromTriplets(R_trip.begin(),R_trip.end(), ass);
		S.setFromTriplets(S_trip.begin(),S_trip.end(), ass);
		Tr.setFromTriplets(Tr_trip.begin(),Tr_trip.end(), ass);
		Ts.setFromTriplets(Ts_trip.begin(),Ts_trip.end(), ass);
		
		this->H=std::move(H);
		
		this->T=std::move(T); this->M=std::move(M); this->N=std::move(N); this->Mp=std::move(Mp); 
		this->Mq=std::move(Mq); this->N=std::move(N); this->P=std::move(P); this->Q=std::move(Q); 
		this->R=std::move(R); this->S=std::move(S); this->Tr=std::move(Tr); this->Ts=std::move(Ts);

		t_material.toc();
		std::cout << "done - " << t_material << " seconds" << std::endl;
	}

	// double add_to_sparse (const double& a, const double& b)
	// {
		// return a+b;
	// }
	
	// double overwrite_to_sparse (const double& a, const double& b)
	// {
		// return b;
	// }
	
	const std::vector<sgnint32_t<int32_t>>& vtf(const int32_t& v_id) const
	{
		if ( vtf_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return vtf_list[uint32_t(abs(v_id))];
	}

	// bool is_conductor(const uint32_t& v)
	// {
		// return vol_material[v] == conductor_id? true : false;
	// }

	const std::vector<sgnint32_t<int32_t>>& fte(const int32_t& f_id) const
	{
		if ( fte_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return fte_list[uint32_t(abs(f_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& ftv(const int32_t& f_id) const
	{
		if ( ftv_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return ftv_list[uint32_t(abs(f_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& etf(const int32_t& e_id) const
	{
		if ( etf_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return etf_list[uint32_t(abs(e_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& etn(const int32_t& e_id) const
	{
		if ( etn_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return etn_list[uint32_t(abs(e_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& nte(const int32_t& n_id) const
	{
		if ( nte_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return nte_list[uint32_t(abs(n_id))];
	}	

	Eigen::Vector3d face_barycenter(const uint32_t& ff)
	{
		if (face_bars.size()!=0)
			return face_bars[ff];
		else
		{
			for (uint32_t f=0; f < surfaces_size(); f++)
			{				
				std::vector<uint32_t> nodes;
				Eigen::Vector3d bc(0,0,0);
				
				for (const auto& signed_ee : fte_list[f])
				{
					uint32_t ee = signed_ee.Val();
					
					for (const auto& signed_nn : etn_list[ee])
					{
						uint32_t nn = signed_nn.Val();
						
						if (!std::binary_search(nodes.begin(),nodes.end(),nn))
						{
							nodes.push_back(nn);
							bc += pts[nn];
						}	
					}
					
				}
				
				face_bars.push_back(bc/3);
			}
			return face_bars[ff];
		}
	}

	Eigen::Vector3d vol_barycenter(const uint32_t& v)
	{
		std::vector<uint32_t> nodes;
		Eigen::Vector3d bc(0,0,0);
		
		for (const auto& signed_ff : vtf_list[v])
		{
			uint32_t ff = signed_ff.Val();
			
			if (ff<surfaces.size())
			{
				for (const auto& signed_ee : fte_list[ff])
				{
					uint32_t ee = signed_ee.Val();
					
					for (const auto& signed_nn : etn_list[ee])
					{
						uint32_t nn = signed_nn.Val();
						
						if (!std::binary_search(nodes.begin(),nodes.end(),nn))
						{
							nodes.push_back(nn);
							bc += pts[nn];
						}	
					}
					
				}
			}
		}
		
		return bc/4;
	}

	Eigen::Vector3d edge_barycenter(const uint32_t& ee)
	{
		if (edge_bars.size()!=0)
			return edge_bars[ee];
		else
		{
			for (uint32_t e = 0; e < edges_size(); e++)
			{
				Eigen::Vector3d bc(0,0,0);
				
				for (const auto& signed_nn : etn_list[e])
				{
					uint32_t nn = signed_nn.Val();
					bc += pts[nn];	
				}
				
				edge_bars.push_back(bc/2);
			}
			
			return edge_bars[ee];
		}
	}

	std::string print_face(const uint32_t& label, const uint32_t& f, bool orient, double r, double g, double b)
	{
		std::ostringstream fr;
		std::set<uint32_t> nodes;
		std::vector<Eigen::Vector3d > n;
		
		for (auto ee : fte(f))
		{
			for (auto nn : etn(abs(ee)))
				nodes.insert(abs(nn));
			
			if (nodes.size() == 3)
				break;
		}
		
		fr << "103.000 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		for (auto nn : nodes)
			n.push_back(pts[nn]);
		
		if (!orient)
			std::swap(n[1],n[2]);
		
		fr << n[0][0] << " " << n[0][1] << " " << n[0][2] << " " ;
		fr << n[1][0] << " " << n[1][1] << " " << n[1][2] << " " ;
		fr << n[2][0] << " " << n[2][1] << " " << n[2][2] << " " ;
			
		fr << std::endl;
		
		return fr.str();
	}

	std::string print_dual_edge(const uint32_t& label, const uint32_t& v, const uint32_t& f, bool orient, double r, double g, double b)
	{
		std::ostringstream fr;
		
		auto fb =  face_barycenter(f);
		auto vb =  vol_barycenter(v);
		
		fr << "102.100 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		if (orient)
			fr << fb[0] << " " << fb[1] << " " << fb[2] << " " << vb[0] << " " << vb[1] << " " << vb[2] << std::endl;
		else
			fr << vb[0] << " " << vb[1] << " " << vb[2] << " " << fb[0] << " " << fb[1] << " " << fb[2] << std::endl;
		
		return fr.str();
	}

	std::vector<std::string> print_dual_face(const uint32_t& label, const uint32_t& e, bool orient, double r, double g, double b)
	{
		// std::ostringstream fr;
		std::vector<std::string> ret;
		auto eb =  edge_barycenter(e);
		double x,y,z;
		x=y=z=0;
		int32_t new_or = orient? 1 : -1;
		
		
		for (auto n : etn_list[e])
		{
			auto segno = n.Sgn();
			x+= new_or*segno*pts[abs(n)][0];
			y+= new_or*segno*pts[abs(n)][1];
			z+= new_or*segno*pts[abs(n)][2];
		}

		Eigen::Vector3d v1 { x, y, z };
		
		for (auto f : etf_list[e])
		{
			auto fb = face_barycenter(abs(f));
			
			for (auto v : ftv_list[abs(f)])
			{
				orient = new_or*f.Sgn()>0;
				auto vb =  vol_barycenter(abs(v));

				Eigen::Vector3d v2 = eb-vb;
				Eigen::Vector3d v3 = fb-vb;
				
				if (v1.dot(v2.cross(v3))>0)
					orient= true;
				else
					orient = false;
					
				
				ret.push_back(print_face(label,orient,vb[0],vb[1],vb[2],eb[0],eb[1],eb[2],fb[0],fb[1],fb[2],r,g,b));
			}
		}
		
		return ret;
	}

	std::string print_face(const uint32_t& f, const int32_t& orient)
	{
		std::ostringstream fr;
		std::set<uint32_t> nodes;
		
		for (auto ee : fte(f))
		{
			for (auto nn : etn(abs(ee)))
				nodes.insert(abs(nn));
			
			if (nodes.size() == 3)
				break;
		}
		
		for (auto nn : nodes)
			fr << nn << " " ;
			
		fr << orient << std::endl;
		
		return fr.str();
	}

	std::string print_edge(const uint32_t& label, const uint32_t& e,
	                       bool orient, double r, double g, double b)
	{
		std::ostringstream fr;
		
		const auto& nn = etn_list[e];
		uint32_t nn_b = abs(*nn.begin());
		uint32_t nn_e = abs(*std::prev(nn.end()));
		Eigen::Vector3d n1 =  pts[nn_b];
		Eigen::Vector3d n2 =  pts[nn_e];
		
		fr << "102.100 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		if (orient)
			fr << n1[0] << " " << n1[1] << " " << n1[2] << " " << n2[0] << " " << n2[1] << " " << n2[2] << std::endl;
		else
			fr << n2[0] << " " << n2[1] << " " << n2[2] << " " << n1[0] << " " << n1[1] << " " << n1[2] << std::endl;
		
		return fr.str();
	}

	std::string print_edge(const uint32_t label, 
	double x1, double y1, double z1, 
	double x2, double y2, double z2)
	{
		std::ostringstream fr;
		
		fr << "102.000 " << label << " " << 0 << " " << 0 << " " << 0 << " 0.0 ";
		fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << std::endl;
		
		return fr.str();
	}

	std::string print_face(const uint32_t label, bool orient, 
	                       double x1, double y1, double z1, 
						   double x2, double y2, double z2, 
						   double x3, double y3, double z3, 
						   double r, double g, double b)
	{
		std::ostringstream fr;
		
		fr << "103.000 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		if (orient)
			fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << " " << x3 << " " << y3 << " " << z3 << std::endl;
		else
			fr << x1 << " " << y1 << " " << z1 << " " << x3 << " " << y3 << " " << z3 << " " << x2 << " " << y2 << " " << z2 << std::endl;
		
		return fr.str();
	}

	uint32_t volumes_size()
	{ 
		return volumes.size(); 
	}
	uint32_t surfaces_size()
	{ 
		return surfaces.size(); 
	}
	uint32_t edges_size()
	{ 
		return edges.size(); 
	}
	uint32_t nodes_size()
	{ 
		return pts.size(); 
	}
	
	private:
	uint32_t									input_line, H_size, Q_size, P_size, N_size, R_size, S_size, M_a_size, M_b_size, M_c_size, U_frac_size, F_frac_size;
	Eigen::SparseMatrix<double> 				C,H,M,Mq,Mp,N,P,Q,R,S,T,Tr,Ts;
	Eigen::VectorXd								U,F,I;
	std::mutex									meshlock;
	std::map<uint32_t,Simulation>				Simulations;
	std::map<uint32_t,Source>					Sources;								/* a std::map works because every time I use the [] operator on an undefined material */
	std::map<uint32_t,Material>					Materials; 								/* (or source), the default constructor makes it empty space (or null source) */
	std::map<uint32_t,BoundaryCondition>		BCs;
	std::map<uint32_t,Mesh>						Meshes;
	std::vector<double>                         CellVolumes;
	std::vector<bool>							is_dirichlet;
	std::vector<uint32_t> 						vol_material, h_index, p_index, n_index, r_index, edge_bids;
	std::vector<uint32_t>				        edge_bcs, face_bcs;
	std::vector<uint8_t>						classify_edges, classify_surfaces;
	std::vector<std::vector<uint32_t>>          associated_volumes;
	std::vector<std::vector<uint32_t>>			edge_src, face_src;
	std::vector<volume_type> 					volumes;
	std::vector<surface_type> 					surfaces;
	std::vector<uint32_t>                       dual_is_fractured, primal_is_fractured;
	std::vector<edge_type> 						edges;
	std::vector<Eigen::Vector3d>	 			pts, edge_bars, face_bars;
	std::vector<cluster_list>    				nte_list, etn_list, etf_list, fte_list, ftv_list, vtf_list;								
	double                                      t_step;
	uint32_t									loaded_mesh_label;
	// std::array<std::vector<uint32_t>,20>		sources_by_label; 						/* Each label has a vector containing all the sources active on that label */
};
