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
 
 //file Solid.hpp
 #include "Utilities.hpp"
 
class Solid
{
	public:
	Solid()
	{
		squareradius = 0;
		material = 0;
		type = "none";
		center = Eigen::Vector3d({0,0,0});
		corner = Eigen::Vector3d({0,0,0});
		size   = Eigen::Vector3d({0,0,0});
	}
	
	//setters
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "type")
		{
			if (std::find(solidtypes.begin(),solidtypes.end(),value) == solidtypes.end())
			{
				std::cout << mesh_throw_preamble;
				MyThrow(input_line,solid_unknown_type);
			}
			type = value;
		}
		else if (param == "material")
			material = std::atoi(value.c_str());
		else if (param == "radius")
		{
			double r = std::stod(value);
			
			if (r <= 0)
			{
				std::cout << mesh_throw_preamble;
				MyThrow(input_line,solid_negative_value);
			}
			
			squareradius = std::pow(r,2);
		}
		else if (param == "center")
		{
			auto  i = value.begin();
			if (*i != '{')
			{
				std::cout << mesh_throw_preamble;
				MyThrow(input_line,coordinates_syntax);
			}
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
					{
						std::cout << mesh_throw_preamble;
						MyThrow(input_line,unbalanced_bracket);
					}
					else 
					{
						if (k < 3)
						{
							center(k)= std::stod(coord);
							k++;
						}
						else
						{
							std::cout << mesh_throw_preamble;
							MyThrow(input_line,too_many_coords);
						}
						
						if (*i == ',')
							i++;
					}
					
				}
				
				if (i == value.end())
				{
					std::cout << mesh_throw_preamble;
					MyThrow(input_line,too_few_coords);
				}
			}
		}
		else if (param == "corner")
		{
			auto  i = value.begin();
			if (*i != '{')
			{
				std::cout << mesh_throw_preamble;
				MyThrow(input_line,coordinates_syntax);
			}
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
					{
						std::cout << mesh_throw_preamble;
						MyThrow(input_line,unbalanced_bracket);
					}
					else 
					{
						if (k < 3)
						{
							corner(k)= std::stod(coord);
							k++;
						}
						else
						{
							std::cout << mesh_throw_preamble;
							MyThrow(input_line,too_many_coords);
						}
						
						if (*i == ',')
							i++;
					}
					
				}
				
				if (i == value.end())
				{
					std::cout << mesh_throw_preamble;
					MyThrow(input_line,too_few_coords);
				}
			}
		}
		else if (param == "size")
		{
			auto  i = value.begin();
			if (*i != '{')
			{
				std::cout << mesh_throw_preamble;
				MyThrow(input_line,coordinates_syntax);
			}
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
					{
						std::cout << mesh_throw_preamble;
						MyThrow(input_line,unbalanced_bracket);
					}
					else 
					{
						if (k < 3)
						{
							size(k)= std::stod(coord);
							if (size(k) <= 0)
							{
								std::cout << mesh_throw_preamble;
								MyThrow(input_line,solid_negative_value);
							}
							k++;
						}
						else
						{
							std::cout << mesh_throw_preamble;
							MyThrow(input_line,too_many_coords);
						}
						
						if (*i == ',')
							i++;
					}
					
				}
				
				if (i == value.end())
				{
					std::cout << mesh_throw_preamble;
					MyThrow(input_line,too_few_coords);
				}
			}
		}
		else
			MyThrow(input_line,solid_unknown_parameter);
	}

	//getters
	bool Rule(const double& x, const double& y, const double& z)
	{
		auto p = Eigen::Vector3d({x,y,z});
		if (type == "sphere")
			return SphereRule(p);
		else if (type == "box")
			return BoxRule(p);
		else if (type == "cylinder")
			return CylRule(p);
		else 
			return false;
	}
	
	bool SphereRule(const Eigen::Vector3d& p)
	{
		Eigen::Vector3d diff = p - center;
		if (diff.dot(diff) <= squareradius)
			return true;
		return false;
	}
	
	bool BoxRule(const Eigen::Vector3d& p)
	{
		auto diff = p - corner;
		if (diff(0) > size(0) || diff(0) < 0)
			return false;
		if (diff(1) > size(1) || diff(1) < 0)
			return false;
		if (diff(2) > size(2) || diff(2) < 0)
			return false;
		return true;
	}
	
	bool CylRule(const Eigen::Vector3d& p)
	{
		if (size.norm() == 0)
			return false;
		
		double dot, dsq;
		Eigen::Vector3d pd = p - center; 
		dot = pd.dot(size);
		
		if( dot < 0.0 || dot > size.norm() )
			return false;

		dsq = pd.squaredNorm() - dot*dot/size.norm();

		if( dsq > squareradius )
			return false;
		
		return true;
	}
	
	const uint32_t& Material(void) const { return material; }
	
	private:
	SolidType type;
	uint32_t material;
	double squareradius;
	Eigen::Vector3d center, corner, size;
};
