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
 
 //file Mesh.hpp
 #include "Utilities.hpp"
 
class Mesh
{
	public:
	Mesh() 
	: file("none"), name("none"), type("none"), mesher("none")
	{
		loaded = false;
		xstep = ystep = zstep = 0;
		xmin = ymin = zmin = 0;
		xmax = ymax = zmax = 0;
		scalefactor = 1;
	}

	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "file")
		{
			std::ifstream test_msh_file(value.c_str());
			if (test_msh_file.is_open())
			{
				test_msh_file.close();
				file = value;
			}
			else
			{
				test_msh_file.close();
				MyThrow(input_line,mesh_inexistent_file);
			}
		}
		else if (param == "name")
			name = value;
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
		else if (param == "xgrid")
		{
			// double xmin,xstep,xmax;
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
						if (k == 0)
						{
							xmin = std::stod(coord);
							k++;
						}
						else if (k == 1)
						{
							xstep = std::stod(coord);
							k++;
						}
						else if (k == 2)
						{
							xmax = std::stod(coord);
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
			
			assert(xmin <= xmax);
			assert(xstep > 0);
		}
		else if (param == "ygrid")
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
						if (k == 0)
						{
							ymin = std::stod(coord);
							k++;
						}
						else if (k == 1)
						{
							ystep = std::stod(coord);
							k++;
						}
						else if (k == 2)
						{
							ymax = std::stod(coord);
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
			
			assert(ymin <= ymax);
			assert(ystep > 0);
		}
		else if (param == "zgrid")
		{
			// double xmin,xstep,xmax;
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
						if (k == 0)
						{
							zmin = std::stod(coord);
							k++;
						}
						else if (k == 1)
						{
							zstep = std::stod(coord);
							k++;
						}
						else if (k == 2)
						{
							zmax = std::stod(coord);
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
			
			assert(zmin <= zmax);
			assert(zstep > 0);
		}
		else if (param == "scalefactor")
			scalefactor = std::stod(value);
		else
			MyThrow(input_line,mesh_unknown_parameter);
	}
	
	const std::string& Name() { return name; }
	const std::string& FileName() { return file; }
	const std::string& GetMeshType() { return type; }
	const std::string& GetMesher()   { return mesher; }
	const double& Scale() { return scalefactor; }
	const double& GetXmin() { return xmin; }
	const double& GetYmin() { return ymin; }
	const double& GetZmin() { return zmin; }
	const double& GetXstep() { return xstep; }
	const double& GetYstep() { return ystep; }
	const double& GetZstep() { return zstep; }
	const double& GetXmax() { return xmax; }
	const double& GetYmax() { return ymax; }
	const double& GetZmax() { return zmax; }
	bool IsLoaded() { return loaded; }
	void Switch() { loaded = !loaded; }
	
	private:
	std::string file, name, type;
	std::string mesher;
	bool loaded;
	double scalefactor;
	double xmin,xmax,ymin,ymax,zmin,zmax,xstep,ystep,zstep,L,volume; //used only when mesh type is cartesian
};

class Refinement
{
	public:
	Refinement()
	{
		xref = 1;
		yref = 1;
		zref = 1;
	}
	
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "vector")
		{
			// double xmin,xstep,xmax;
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
						if (k == 0)
						{
							xref = std::stod(coord);
							k++;
						}
						else if (k == 1)
						{
							yref = std::stod(coord);
							k++;
						}
						else if (k == 2)
						{
							zref = std::stod(coord);
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

			assert(xref >= 1);
			assert(yref >= 1);
			assert(zref >= 1);
		}
		else if (param == "x")
		{
			xref = std::stod(value);
			assert(xref >= 1);
		}
		else if (param == "y")
		{
			yref = std::stod(value);
			assert(yref >= 1);
		}
		else if (param == "z")
		{
			zref = std::stod(value);
			assert(zref >= 1);
		}
		else
			MyThrow(input_line,ref_unknown_parameter);
	}
	
	const double& X(void) { return xref; }
	const double& Y(void) { return yref; }
	const double& Z(void) { return zref; }
	
	private:
	double xref,yref,zref;
};
