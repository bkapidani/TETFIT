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
#include <set>
#include <limits>
#include <sstream>
#include <future>
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

namespace parser {

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
} //namespace priv

using cluster_list     		= std::vector<sgnint32_t<int32_t>>; 
using volume_type 			= std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>;
using surface_type 			= std::tuple<uint32_t,uint32_t,uint32_t>;
using label_surface_type 	= std::pair<surface_type,uint32_t>;
using edge_type 			= std::tuple<uint32_t,uint32_t>;
using label_edge_type 		= std::pair<edge_type,uint32_t>;
using label_node_type 		= std::pair<uint32_t,uint32_t>;
using tm_tuple 				= std::tuple< volume_type, uint32_t >;
using sm_tuple 				= std::tuple< surface_type, uint32_t >;
using double_triplet		= Eigen::Triplet<double,uint32_t>;

typedef std::string				Primitive;
typedef std::string				Sourcetype;
typedef std::string				Direction;
typedef std::string				BaseFunction;
typedef std::string				BoundaryConditionType;
typedef std::string				Profile;
typedef std::string				Meshtype;
typedef std::string				OutputMode;
typedef double					Amplitude;
typedef double 					Frequency;
typedef double 					WaveNumber;
typedef double					Duration;
typedef std::array<double,4>	SpaceTimePoint;
typedef std::array<double,3> 	WaveVector;


//lists of allowed string constants
const std::vector<Primitive>							definables		= {"material","source","mesh","bc","simulation","geometry"};
const std::vector<Sourcetype>   						sourcetypes   	= { "sine", "gaussian", "dc" };
const std::vector<Profile>   							profiles   		= { "point", "square", "cube", "domain"};
const std::vector<Direction>    						directions    	= { "x", "y", "z" };
const std::vector<BaseFunction> 						basefunctions 	= { "sin", "cos" };
const std::vector<BoundaryConditionType>				bctypes			= { "pec", "pmc", "pml" };
const std::vector<Meshtype>								meshtypes		= { "tetrahedral", "cartesian", "none"};
const std::vector<OutputMode>							outputmodes		= { "silo", "probe"};

template<typename T>
void sort_unique(std::vector<T>& v) //useful as stand-alone
{
	std::sort(v.begin(), v.end());
	auto uniq_iter = std::unique(v.begin(), v.end());
	v.erase(uniq_iter, v.end());
}

class BoundaryCondition
{
	public:
	BoundaryCondition()
	: type("pec")
	{
		thickness=0;
	}
	
	void SetParam(std::string param, std::string value)
	{
		if (param == "type")
		{
			// Temporary error message
			if (value == "pml")
				throw std::runtime_error(std::string("Sorry, PML not implemented yet, getting there!"));
			else if (value == "pmc")
				throw std::runtime_error(std::string("Sorry, PMC not implemented yet, getting there!"));
			else if (value == "pec")
			{
				type = value;
				this->val=0;
			}
			else
				throw std::runtime_error(std::string("Input file error: Unrecognized boundary condition type! Available: pec, pmc, pml"));				
		}
		else if (param == "thickness")
			thickness = std::stod(value);
		else	
			throw std::runtime_error(std::string("Input file error: Unrecognized boundary condition parameter! Available: type"));	
	}
	const std::string& 		GetBCType(void) { return type; }
	const double& 			GetThickness(void) { return thickness; }
	const double&			GetValue(void) { return val; }
	
	private:
	std::string type;
	double thickness, val;
};

class Source
{
	public:
	Source()
	{
		st    	= "dc";
		dir   	= "x";
		amp   	= 0;
		freq  	= 0;
		kvec[0] = kvec[1] = kvec[2] = 0;
		prof    = "point";
		bfuncs[0] = bfuncs[1] = bfuncs[2] = "cos";
		center_coords = {0,0,0};
		surface_label   = 0;
	}
	
	void SetParam(std::string param, std::string value)
	{
		if (param == "direction")
		{
			if (std::find(directions.begin(),directions.end(),value) == directions.end())
				throw std::runtime_error(std::string("Input file error: Unrecognized direction! Available: x, y, z"));
			else
				this->dir = value;
		}
		else if (param == "type")
		{
			if (std::find(sourcetypes.begin(),sourcetypes.end(),value) == sourcetypes.end())
				throw std::runtime_error(std::string("Input file error: Unrecognized source type! Available: Gaussian, sinewave, dc"));
			else
				this->st = value;
		}
		else if (param == "profile")
		{
			if (std::find(profiles.begin(),profiles.end(),value) == profiles.end())
				throw std::runtime_error(std::string("Input file error: Unrecognized source profile! Available: point, square, cube, domain"));
			else
				this->prof = value;
		}
		else if (param == "center")
		{
			auto  i = value.begin();
			if (*i != '{')
				throw std::runtime_error(std::string("Input file error: center coordinates must be inside braces {...,..,..}"));
			else
			{
				uint8_t k=0;
				i++;
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					while (*i != ',' && *i != '}' && i != value.end());
					{
						coord.push_back(*i);
						i++;
					}
					if (i == value.end())
						throw std::runtime_error(std::string("Input file error: unbalanced bracket"));
					else 
					{
						if (k < 3)
						{
							center_coords[k]= std::stod(coord);
							k++;
						}
						else
							throw std::runtime_error(std::string("Input file error: Maximum of three coordinates!"));
						
						if (*i = ',')
							i++;
					}
				}
				
				if (i == value.end())
					throw std::runtime_error(std::string("Input file error: Undefined end to list of coordinates!"));
			}
		}
		else if (param == "wavevector")
		{
			auto  i = value.begin();
			if (*i != '{')
				throw std::runtime_error(std::string("Input file error: wavevector coordinates must be inside braces {...,..,..}"));
			else
			{
				uint8_t k=0;
				i++;
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					while (*i != ',' && *i != '}' && i != value.end());
					{
						coord.push_back(*i);
						i++;
					}
					if (i == value.end())
						throw std::runtime_error(std::string("Input file error: unbalanced bracket"));
					else 
					{
						switch(k)
						{
							case 0 : 
							{
								kvec[0] = std::stod(coord);
								k++;
								break;
							}
							case 1 : 
							{
								kvec[1] = std::stod(coord);
								k++;
								break;
							}
							case 2 : 
							{
								kvec[2] = std::stod(coord);
								k++;
								break;
							}
							default : 
							{	
								throw std::runtime_error(std::string("Input file error: Maximum of three coordinates!"));
								break;
							}
						}
						
						if (*i = ',')
							i++;
					}
					
				}
				
				if (i == value.end())
					throw std::runtime_error(std::string("Input file error: Undefined end to list of coordinates!"));
			}
		}
		else if (param == "basefunctionvector")
		{
			auto  i = value.begin();
			if (*i != '{')
				throw std::runtime_error(std::string("Input file error: base-function vector coordinates must be inside braces {...,..,..}"));
			else
			{
				uint8_t k=0;
				i++;
				while (*i != ',' && *i != '}' && i != value.end())
				{
					std::string coord;
					while (*i != ',' && *i != '}' && i != value.end());
					{
						if (*i != ' ')	
							coord.push_back(*i);
						i++;
					}
					if (i == value.end())
						throw std::runtime_error(std::string("Input file error: unbalanced bracket"));
					else 
					{
						switch(k)
						{
							case 0 : 
							{
								bfuncs[0] = coord;
								k++;
								break;
							}
							case 1 : 
							{
								bfuncs[1] = coord;
								k++;
								break;
							}
							case 2 : 
							{
								bfuncs[2] = coord;
								k++;
								break;
							}
							default : 
							{	
								throw std::runtime_error(std::string("Input file error: Maximum of three coordinates!"));
								break;
							}
						}
						
						if (*i = ',')
							i++;
					}
				}
				
				if (i == value.end())
					throw std::runtime_error(std::string("Input file error: Undefined end to list of coordinates!"));
			}
		}
		else if (param == "Xbasefunction")
		{
			if (std::find(basefunctions.begin(),basefunctions.end(),value) == basefunctions.end())
				throw std::runtime_error(std::string("Input file error: Unrecognized base function! Available: sin, cos"));
			else
				this->bfuncs[0] = value;
		}
		else if (param == "Ybasefunction")
		{
			if (std::find(basefunctions.begin(),basefunctions.end(),value) == basefunctions.end())
				throw std::runtime_error(std::string("Input file error: Unrecognized base function! Available: sin, cos"));
			else
				this->bfuncs[1] = value;
		}
		else if (param == "Zbasefunction")
		{
			if (std::find(basefunctions.begin(),basefunctions.end(),value) == basefunctions.end())
				throw std::runtime_error(std::string("Input file error: Unrecognized base function! Available: sin, cos"));
			else
				this->bfuncs[2] = value;
		}
		else if (param == "amplitude")
			this->amp = std::stod(value);
		else if (param == "frequency")
			this->freq = std::stod(value);
		else if (param == "kx")
			this->kvec[0] = std::stod(value);
		else if (param == "ky")
			this->kvec[1] = std::stod(value);
		else if (param == "kz")
			this->kvec[2] = std::stod(value);
		else if (param == "surface")
			this->surface_label = std::stod(value);
		else
			throw std::runtime_error(std::string("Input file error: unrecognized parameter for source!"));
		
	}
	
	const Sourcetype& 	GetSourceType(void) { return st; } 
	const Direction& 	GetDirection(void) { return dir; }
	const WaveNumber&	Getkx(void) { return kvec[0]; }
	const WaveNumber&	Getky(void) { return kvec[1]; }
	const WaveNumber&	Getkz(void) { return kvec[2]; }
	const Frequency&	GetFreq(void) { return freq; }
	const Amplitude&	GetAmp(void) { return amp; }
	
	Eigen::Vector3d Compute(SpaceTimePoint p)
	{
		double ret = amp*cos(2*PI*freq*p[3]);
		
		for (uint8_t j=0; j<3; j++)
		{
			if (bfuncs[j] == "sin")
				ret *= sin(2*PI*kvec[j]*p[j]);
			else
				ret *= cos(2*PI*kvec[j]*p[j]);
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
	Direction dir;
	WaveVector kvec;
	std::array<BaseFunction,3> bfuncs;
	std::vector<double> center_coords;
	uint32_t surface_label;
	// double (*foo)();
};

class Material
{
	public:
	Material()
	{
		//All values initialized to vacuum parameters
	}
	
	//setters
	void SetParam(std::string param, std::string value)
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
			throw std::runtime_error(std::string("Input file error: Unrecognized material parameter! Available: epsilon, mu, sigma, chi"));
	}

	//getters
	double Epsilon(void) { return epsilon; }
	double Mu(void) { return mu; }
	double Sigma(void) { return sigma; }
	double Chi(void) { return chi; }
	
	private:
	double epsilon = epsilon0;
	double sigma = 0;
	double mu = mu0;
	double chi = 0;
};

class Mesh
{
	public:
	Mesh() 
	: file("cube.mesh"), type("none"), mesher("netgen")
	{}

	void SetParam(std::string param, std::string value)
	{
		if (param == "file")
			file = value;
		else if (param == "type")
		{
			if (std::find(meshtypes.begin(),meshtypes.end(),value) == meshtypes.end())
				throw std::runtime_error(std::string("Input file error: undefined mesh type! Available: tetrahedral, cartesian "));
			type = value;
		}
		else if (param == "mesher")
			mesher = value;
		else if (param == "xstep")
			xstep = std::stod(value);
		else if (param == "ystep")
			ystep = std::stod(value);
		else if (param == "zstep")
			zstep = std::stod(value);
		else
			throw std::runtime_error(std::string("Input file error: undefined mesh parameter! Available: file, type"));
	}
	
	const std::string& GetFileName() { return file; }
	const std::string& GetMeshType() { return type; }
	const std::string& GetMesher()   { return mesher; }
	
	private:
	std::string file;
	std::string type;
	std::string mesher;
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
	
	void SetParam(std::string param, std::string value)
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
				throw std::runtime_error(std::string("Input file error: undefined output mode type! Available: silo, probe"));
			mode = value;
		}
		else
			throw std::runtime_error(std::string("Input file error: undefined simulation parameter! Available: source, mesh, duration"));
	}
	
	private:
	Duration d;
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
		std::ifstream readFile(inputfile.c_str());
		
		bool in_definition = false;
		bool mesh_added = false;
		
		
		char action[10], token[20], value[64];
		std::string thing_being_defined;
		uint32_t definition_label;
		
		std::cout << "ciao!" << std::endl;
		
		while(getline(readFile,line))
		{
			std::cout << "ciao!" << std::endl;
			auto c_line = line.c_str();
			if (line.size()>0    &&
				c_line[0] != '#' &&       /* lines that begin with '#' are comments         */
				c_line[0] != '\n'  )      /* empty lines for better readability are allowed */
			{
				sscanf(c_line,"%s %s %s",action,token,value);
				std::string instr(action), tok(token), val(value);
				
				if (!in_definition)
				{
					if (instr != "DEFINE")
						throw std::runtime_error(std::string("Input file error: define something before setting variables"));
					else 
					{
						if (std::find(definables.begin(),definables.end(),tok) == definables.end())
							throw std::runtime_error(std::string("Input file error: can only define material, mesh, boundary condition or source"));
						else
						{
							thing_being_defined = tok;
							definition_label = std::stod(val);
							in_definition = true;
						}
					}
				}
				else if (instr == "END")
				{
					if (tok != thing_being_defined || std::atoi(value) != definition_label)
						throw std::runtime_error(std::string("Input file error: ending non defined definition"));
					else 
						in_definition = false;
				}
				else if (instr != "SET")
					throw std::runtime_error(std::string("Input file error: Unknown instruction inside define block"));
				else
				{
					if (thing_being_defined == "material")
						Materials[definition_label].SetParam(tok,val);
					else if (thing_being_defined == "bc")
						BCs[definition_label].SetParam(tok,val);
					else if (thing_being_defined == "mesh")
						Meshes[definition_label].SetParam(tok,val);
					else if (thing_being_defined == "source")
						Sources[definition_label].SetParam(tok,val);
					else if (thing_being_defined == "simulation")
						Simulations[definition_label].SetParam(tok,val);
					else
						throw std::runtime_error(std::string("Input file error: can only define material, mesh, boundary condition or source"));
				}
			}
		}
		
		if (in_definition)
			throw std::runtime_error(std::string("Input file error: File ended unexpectedly before ending definition"));
		if (!mesh_added)
			throw std::runtime_error(std::string("Input file error: No mesh file was added!"));

		timecounter t_read;
		t_read.tic();
		ReadMesh((*Meshes.begin()).second.GetFileName());
		t_read.toc();
		std::cout << "Loading complex took: " << t_read << " s" << std::endl;
		
		readFile.close();
	}
	
	/*Discretization(std::string meshfile)
	{
		Epsilon[1]=epsilon0;
		Mu[1]=mu0;
		Sigma[1]=0.02;
		Chi[1]=0;
		freq=5e9;
		timecounter t_read;
		t_read.tic();
		ReadMesh(meshfile);
		t_read.toc();
		std::cout << "Loading complex took: " << t_read << " s" << std::endl;
	}*/

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

	bool ReadMesh(const std::string& _filename)
	{	
		timecounter tc, tctot;
		
		/* Open file */
		if (_filename.size() == 0)
		{
			std::cout << "Invalid mesh file name" << std::endl;
			return false;
		}
		
		uint32_t	lines, linecount;
		
		mapped_file mf(_filename);
		
		std::cout << " * * * Reading NETGEN format mesh * * * " << std::endl;
		
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
			
			
			Eigen::Vector3d point = { std::get<0>(t),std::get<1>(t),std::get<2>(t) };		
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
		
		// vol_signs.reserve(lines);
		// volumes.reserve(lines);
		// vol_material.reserve(lines);
		
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

			auto vol_vol = v1.dot(v2.cross(v3))/double(6)>0;
			int32_t sgn  = vol_vol? 1 : -1;
			vol_signs.push_back(sgn);
			CellVolumes.push_back(vol_vol);
			vol_material.push_back(std::get<1>(tet));
		}
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
			
			tripletList.push_back(double_triplet(k,double(abs(e1)), 1));
			tripletList.push_back(double_triplet(k,double(abs(e2)),-1));
			tripletList.push_back(double_triplet(k,double(abs(e3)), 1));
			
			auto vols = ftv_list[k];
			
			bool recombine = true;
			switch (vols.size()) 
			{
				case 2: 
				{
					auto vol1= abs(vols[0]);
					auto vol2= abs(vols[1]);
					
					if (Materials[vol_material[vol1]].Chi() != 0 || Materials[vol_material[vol2]].Chi() != 0)
					{
						//We have magnetic losses
						if (Materials[vol_material[vol1]].Chi() != Materials[vol_material[vol2]].Chi())
							recombine=false;
						else if (Materials[vol_material[vol1]].Mu() != Materials[vol_material[vol2]].Mu())
							recombine=false;
					}
					
					break;
				}
				case 1:
				{
					/*Don't know what to make of this, yet*/
					break;
				}
				case 0:
				{
					throw std::invalid_argument("Conductor boundary cannot be on mesh boundary!");
					break;
				}
			}

			recombine_faces.push_back(recombine);
		}
		
		std::vector<uint32_t>().swap(e_labels);
		nte_list.resize(pts.size());
		
		C.setFromTriplets(tripletList.begin(),tripletList.end());
		this->C=std::move(C);
		
		// etn_list.reserve(edges.size());
		// physical_nodes.resize(pts.size(),0);
		
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
					while (itor2 != associated_volumes[p0].end())
					{

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
			
			recombine_edges.push_back(recombine);
				
		}
		
		tc.toc();
		
		std::cout << "done - " << tc << " seconds" << std::endl;
		
		/************************ Read boundary surfaces ************************/
		linecount = 0;
		auto num_of_tets=lines;
		lines = strtot<uint32_t>(endptr, &endptr);
		
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
				std::cout << "Reading triangle: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}
			
			auto t = parser::read_triangle_line<uint32_t>(endptr, &endptr);

			
			uint32_t       p0( std::get<1>(t) );
			uint32_t       p1( std::get<2>(t) );
			uint32_t       p2( std::get<3>(t) );
			uint32_t       bid( std::get<0>(t) );
			
			surface_type   tri( p0, p1, p2 );
			
			// if (!physical_surfaces[bid].size())
				// physical_surfaces[bid].resize(surfaces.size(),0);

			// auto itor = std::lower_bound(surfaces.begin(),surfaces.end(),tri);
			// physical_surfaces[bid][std::distance(surfaces.begin(),itor)]++;
			
			linecount++;
		}
		tc.toc();
		
		std::cout << "Reading triangle: " << linecount;
		std::cout << "/" << lines  << " - " << tc << " seconds"  << std::endl;

		tctot.toc();
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
        for (auto itor = 0; itor < volumes_size(); itor++)
		{
			auto vol = CellVolumes[itor];
			auto vol_domain = vol_material[itor];
			double eps_vol = Materials[vol_domain].Epsilon();
			double mu_vol = Materials[vol_domain].Mu();
			double c = 1/sqrt(eps_vol*mu_vol);
			auto fareas = primal_area_vectors(itor);
			
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
		
		std::cout << "Minimum diameter: " << max_h << std::endl;
		
		return ret;
	}
	
	Eigen::VectorXd power_method_iteration(Eigen::VectorXd b )
	{
		b = b*(1/b.norm());
		auto cb = C*b;
		auto nicb = Ni*cb;
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
	
	bool ConstructMaterialMatrices(void)
	{
		Eigen::MatrixXd local_E = Eigen::MatrixXd::Zero(50,50);
		Eigen::MatrixXd local_S = Eigen::MatrixXd::Zero(50,50);
		t_step = estimate_time_step_bound();
		
		std::vector<bool> mu_computed(volumes_size(),false);
		
		for (uint32_t nid=0; nid< pts.size(); nid++ )
		{
			// timecounter t_find;
			// t_find.tic();
			std::set<uint32_t> global_i;
			
			// auto nid= _mesh->find_id(*nn);				
			auto vols = associated_volumes[nid];
			
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

				if (!mu_computed[vv])
				{
					double mu_vol = Materials[vol_material[vv]].Mu();
					auto face_vecs = dual_area_vectors(vv);
					auto fids = vtf_list[vv];
					double coeff=36.0/(fabs(CellVolumes[vv])*(1/mu_vol));
					mu_computed[vv]=true;
					
					Eigen::Matrix4d mu;

					mu(0,0)=((face_vecs[0]).dot(face_vecs[0])+(face_vecs[2]).dot(face_vecs[2])+(face_vecs[3]).dot(face_vecs[3]))*coeff;
					mu(1,1)=((face_vecs[0]).dot(face_vecs[0])+(face_vecs[1]).dot(face_vecs[1])+(face_vecs[4]).dot(face_vecs[4]))*coeff;
					mu(2,2)=((face_vecs[1]).dot(face_vecs[1])+(face_vecs[2]).dot(face_vecs[2])+(face_vecs[5]).dot(face_vecs[5]))*coeff;
					mu(3,3)=((face_vecs[3]).dot(face_vecs[3])+(face_vecs[4]).dot(face_vecs[4])+(face_vecs[5]).dot(face_vecs[5]))*coeff;

					mu(0,1)=(-(face_vecs[1]).dot(face_vecs[2])-(face_vecs[3]).dot(face_vecs[4]))*coeff;
					mu(0,2)=(-(face_vecs[0]).dot(face_vecs[1])+(face_vecs[3]).dot(face_vecs[5]))*coeff;
					mu(0,3)=( (face_vecs[0]).dot(face_vecs[4])+(face_vecs[2]).dot(face_vecs[5]))*coeff;

					mu(1,0)=mu(0,1);
					mu(1,2)=(-(face_vecs[0]).dot(face_vecs[2])-(face_vecs[4]).dot(face_vecs[5]))*coeff;
					mu(1,3)=( (face_vecs[0]).dot(face_vecs[3])-(face_vecs[1]).dot(face_vecs[5]))*coeff;
					
					mu(2,0)=mu(0,2);
					mu(2,1)=mu(1,2);
					mu(2,3)=(-(face_vecs[2]).dot(face_vecs[3])-(face_vecs[1]).dot(face_vecs[4]))*coeff;
					
					mu(3,0)=mu(0,3);
					mu(3,1)=mu(1,3);
					mu(3,2)=mu(2,3);

					/*Old version*/
					
					// mu(0,0)=(dot(face_vecs[0],face_vecs[0])+dot(face_vecs[2],face_vecs[2])+dot(face_vecs[3],face_vecs[3]))*coeff;
					// mu(1,1)=(dot(face_vecs[0],face_vecs[0])+dot(face_vecs[1],face_vecs[1])+dot(face_vecs[4],face_vecs[4]))*coeff;
					// mu(2,2)=(dot(face_vecs[1],face_vecs[1])+dot(face_vecs[2],face_vecs[2])+dot(face_vecs[5],face_vecs[5]))*coeff;
					// mu(3,3)=(dot(face_vecs[3],face_vecs[3])+dot(face_vecs[4],face_vecs[4])+dot(face_vecs[5],face_vecs[5]))*coeff;

					// mu(0,1)=(-dot(face_vecs[1],face_vecs[2])-dot(face_vecs[3],face_vecs[4]))*coeff;
					// mu(0,2)=(-dot(face_vecs[0],face_vecs[1])+dot(face_vecs[3],face_vecs[5]))*coeff;
					// mu(0,3)=( dot(face_vecs[0],face_vecs[4])+dot(face_vecs[2],face_vecs[5]))*coeff;

					// mu(1,0)=mu(0,1);
					// mu(1,2)=(-dot(face_vecs[0],face_vecs[2])-dot(face_vecs[4],face_vecs[5]))*coeff;
					// mu(1,3)=( dot(face_vecs[0],face_vecs[3])-dot(face_vecs[1],face_vecs[5]))*coeff;
					
					// mu(2,0)=mu(0,2);
					// mu(2,1)=mu(1,2);
					// mu(2,3)=(-dot(face_vecs[2],face_vecs[3])-dot(face_vecs[1],face_vecs[4]))*coeff;
					
					// mu(3,0)=mu(0,3);
					// mu(3,1)=mu(1,3);
					// mu(3,2)=mu(2,3);
					
					auto local_Ni = mu.inverse();
					
					if (Materials[vol_material[vv]].Chi() != 0)
					{
						//...
					}
					
				}
				
				uint32_t k=0;
				for (auto vn : vol_nodes)
				{
					if (vn == nid)
						break;
					else
						k++;
				}
				
				
				if (k==0)
				{
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[0]);
					it_gi=ret.first;
					
					it_gi = global_i.insert(it_gi,edgs_l2g[2]);

					it_gi = global_i.insert(it_gi,edgs_l2g[3]);
					
				}
				else if (k==1)
				{
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[0]);
					it_gi=ret.first;
					
					it_gi = global_i.insert(it_gi,edgs_l2g[1]);

					it_gi = global_i.insert(it_gi,edgs_l2g[4]);			
				}
				else if (k==2)
				{
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[1]);
					it_gi=ret.first;
					
					it_gi = global_i.insert(it_gi,edgs_l2g[2]);

					it_gi = global_i.insert(it_gi,edgs_l2g[5]);
				}
				else if (k==3)
				{
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[3]);
					it_gi=ret.first;
					it_gi = global_i.insert(it_gi,edgs_l2g[4]);
					it_gi = global_i.insert(it_gi,edgs_l2g[5]);			
					
				}
				
			}
			
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
				double sigma_vol = Materials[vol_domain].Sigma();
				
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
				
				
				if (k==0)
				{
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[0]);
					it_gi=ret.first;							
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = global_i.insert(it_gi,edgs_l2g[2]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = global_i.insert(it_gi,edgs_l2g[3]);						
					e3 = std::distance(global_i.begin(),it_gi);					
					
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
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[0]);
					it_gi=ret.first;							
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = global_i.insert(it_gi,edgs_l2g[1]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = global_i.insert(it_gi,edgs_l2g[4]);
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
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[1]);
					it_gi=ret.first;							
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = global_i.insert(it_gi,edgs_l2g[2]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = global_i.insert(it_gi,edgs_l2g[5]);							
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
					std::set<uint32_t>::iterator it_gi;
					
					auto ret = global_i.insert(edgs_l2g[3]);
					it_gi=ret.first;							
					e1 = std::distance(global_i.begin(),it_gi);
					
					it_gi = global_i.insert(it_gi,edgs_l2g[4]);
					e2 = std::distance(global_i.begin(),it_gi);

					it_gi = global_i.insert(it_gi,edgs_l2g[5]);							
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
				
			}
			
			auto local_E_size = global_i.size();
			local_E.resize(local_E_size,local_E_size);
			local_S.resize(local_E_size,local_E_size);
			
			if (sigma_node != 0)
				local_S*=0.5*t_step;
			
			uint32_t jj,kk;
			jj=kk=0;
			
			/*This goes on its own if I don't want to already put in boundary conditions in the M_eps matrix*/
			/*for (auto j = global_i.begin(); j != global_i.end(); j++)
			{
				std::tuple<uint32_t,uint32_t,T> ct, cs;
				
				if (timestep_version != 3)
				{
					ct = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*j),ctt_size+jj,1);
					ctt_entries.push_back(ct);
				}
				
				if (is_dirichlet.at(*j))
				{
					local_E.row(jj).fill(0);
					local_E.col(jj).fill(0);
					local_E(jj,jj)=1;
					
					if (sigma_node != 0)
					{
						local_S.row(jj).fill(0);
						local_S.col(jj).fill(0);
					}
					
					if (timestep_version == 3)
						ct = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*j),uint32_t(*j),1);
					else
						ct = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+jj,ctt_size+jj,1);
					
					
					Et_entries.push_back(ct);
				}
				else
				{
					for (auto k = j; k != global_i.end(); k++)
					{
						if (is_dirichlet.at(*k) )
							;
						else if (local_E(jj,kk)!=0 && !sigma_node)
						{
							if (timestep_version == 3)
								ct = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*j),uint32_t(*k),double(local_E(jj,kk)));
							else
								ct = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+jj,ctt_size+kk,double(local_E(jj,kk)));

							Et_entries.push_back(ct);
							
							if (jj!=kk)
							{
								if (timestep_version == 3)
									ct = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*k),uint32_t(*j),double(local_E(jj,kk)));
								else
									ct = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+kk,ctt_size+jj,double(local_E(jj,kk)));
								
								Et_entries.push_back(ct);
							}
						}
						else if (sigma_node &&  (local_E(jj,kk)!=0 || local_S(jj,kk)!=0) )
						{
							
							if (timestep_version == 3)
								cs = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*j),uint32_t(*k),double(local_E(jj,kk)-local_S(jj,kk)));
							else
								cs = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+jj,ctt_size+kk,double(local_E(jj,kk)-local_S(jj,kk)));

							Et_entries.push_back(cs);

							if (jj!=kk)
							{
								if (timestep_version == 3)
									cs = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*k),uint32_t(*j),double(local_E(jj,kk)-local_S(jj,kk)));
								else
									cs = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+kk,ctt_size+jj,double(local_E(jj,kk)-local_S(jj,kk)));
								
								Et_entries.push_back(cs);
							}
							

						}
						kk++;
					}
				}
				jj++;
				kk=jj;
			}*/
			
			// t_find.toc();			
			
			// t_find.tic();				
			auto local_Id = Eigen::MatrixXd::Identity(local_E_size,local_E_size);
			auto local_H = (local_E+local_S).llt().solve(local_Id);
			// auto local_H = (local_E+local_S).inverse();
			
			
			// t_find.tic();

			/*jj=kk=0;
				
			for (auto j = global_i.begin(); j != global_i.end(); j++)
			{
				for (auto k = j; k != global_i.end(); k++)
				{
					if (local_H(jj,kk)!=0)
					{
						std::tuple<uint32_t, uint32_t, double> ct;
						
						if (timestep_version == 3)
							ct = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*j),uint32_t(*k),double(local_H(jj,kk)));
						else
							ct = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+jj,ctt_size+kk,double(local_H(jj,kk)));
						
						Ht_entries.push_back(ct);
						
						if (jj != kk)
						{
							if (timestep_version == 3)
								ct = std::make_tuple<uint32_t, uint32_t, double>(uint32_t(*k),uint32_t(*j),double(local_H(jj,kk)));
							else
								ct = std::make_tuple<uint32_t, uint32_t, double>(ctt_size+kk,ctt_size+jj,double(local_H(jj,kk)));
							
							Ht_entries.push_back(ct);
						}
					}
					kk++;
				}
				jj++;
				kk=jj;
			}

			
			ctt_size+=local_E_size;*/
			
			local_E = Eigen::MatrixXd::Zero(50,50);
			local_S = Eigen::MatrixXd::Zero(50,50);
			
			// t_find.toc();
		}
		
	}

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

	Eigen::Vector3d face_barycenter(const uint32_t& f)
	{
		std::vector<uint32_t> nodes;
		Eigen::Vector3d bc(0);
		
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
		
		return bc/3;
	}

	Eigen::Vector3d vol_barycenter(const uint32_t& v)
	{
		std::vector<uint32_t> nodes;
		Eigen::Vector3d bc(0);
		
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

	Eigen::Vector3d edge_barycenter(const uint32_t& e)
	{
		Eigen::Vector3d bc(0);
		
		for (const auto& signed_nn : etn_list[e])
		{
			uint32_t nn = signed_nn.Val();
			bc += pts[nn];	
		}
		
		return bc/2;
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

	std::string print_edge(const uint32_t& label, const uint32_t& e, bool orient, double r, double g, double b)
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

	std::string print_edge(const uint32_t label, double x1, double y1, double z1, double x2, double y2, double z2)
	{
		std::ostringstream fr;
		
		fr << "102.000 " << label << " " << 0 << " " << 0 << " " << 0 << " 0.0 ";
		fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << std::endl;
		
		return fr.str();
	}

	std::string print_face(const uint32_t label, bool orient, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double r, double g, double b)
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
	
	// private:
	Eigen::SparseMatrix<double> 				C,H,E,Ni;
	std::map<uint32_t,Simulation>				Simulations;
	std::map<uint32_t,Source>					Sources;							// a std::map works because every time I use the [] operator on an undefined material
	std::map<uint32_t,Material>					Materials; 							// (or source), the default constructor makes it empty space (or null source)
	std::map<uint32_t,BoundaryCondition>		BCs;
	std::map<uint32_t,Mesh>						Meshes;
	// std::map<uint32_t,double>                   Epsilon,Mu,Sigma,Chi;			
	std::vector<double>                         CellVolumes;
	std::vector<uint32_t> 						vol_material;
	std::vector<bool>							recombine_edges, recombine_faces;
	std::vector<std::vector<uint32_t>>          associated_volumes;
	std::vector<volume_type> 					volumes;
	std::vector<surface_type> 					surfaces;
	std::vector<edge_type> 						edges;
	std::vector<Eigen::Vector3d>	 			pts;
	/* node -> cluster of edge IDs around it */
	std::vector<cluster_list>    				nte_list;
	std::vector<cluster_list>    				etn_list;	
	/* edge -> cluster of face IDs around it */
	std::vector<cluster_list>        			etf_list;
	std::vector<cluster_list>        			fte_list;	
	/* triangle -> cluster of volume IDs around it */
	std::vector<cluster_list>  					ftv_list;
	std::vector<cluster_list> 					vtf_list;
	double                                      freq, t_step;
};

