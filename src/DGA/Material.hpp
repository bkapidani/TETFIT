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
 
 //file Material.hpp
 #include "Utilities.hpp"
 
class Material
{
	public:
	Material()
	{
		epsilon_scalar = epsilon0;
		sigma_scalar   = 0;
		mu_scalar      = mu0;
		chi_scalar     = 0;
		epsilon = epsilon0*Eigen::Matrix3d::Identity();
		sigma = 0*Eigen::Matrix3d::Identity();
		mu = mu0*Eigen::Matrix3d::Identity();
		chi = 0*Eigen::Matrix3d::Identity();
		pml = true;
	}
	
	//setters
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		pml = false;
		if (param == "epsilon")
		{
			epsilon_scalar = epsilon0*std::stod(value);
			epsilon        = epsilon_scalar*Eigen::Matrix3d::Identity();
		}
		else if (param == "epsilon_xy")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			epsilon(0,1) = epsilon(1,0) = std::stod(value);
		}
		else if (param == "epsilon_xz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			epsilon(0,2) = epsilon(2,0) = std::stod(value);
		}
		else if (param == "epsilon_yz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			epsilon(1,2) = epsilon(2,1) = std::stod(value);
		}
		else if (param == "mu")
		{
			mu_scalar = mu0*std::stod(value);
			mu		  = mu_scalar*Eigen::Matrix3d::Identity();
		}
		else if (param == "mu_xy")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			mu(0,1) = mu(1,0) = std::stod(value);
		}
		else if (param == "mu_xz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			mu(0,2) = mu(2,0) = std::stod(value);
		}
		else if (param == "mu_yz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			mu(1,2) = mu(2,1) = std::stod(value);
		}
		else if (param == "sigma")
		{
			sigma_scalar = std::stod(value);
			sigma 		 = sigma_scalar*Eigen::Matrix3d::Identity();
		}
		else if (param == "sigma_xy")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			sigma(0,1) = sigma(1,0) = std::stod(value);
		}
		else if (param == "sigma_xz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			sigma(0,2) = sigma(2,0) = std::stod(value);
		}
		else if (param == "sigma_yz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			sigma(1,2) = sigma(2,1) = std::stod(value);
		}
		else if (param == "chi")
		{
			chi_scalar = std::stod(value);
			chi 	   = chi_scalar*Eigen::Matrix3d::Identity();
		}
		else if (param == "chi_xy")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			chi(0,1) = chi(1,0) = std::stod(value);
		}
		else if (param == "chi_xz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			chi(0,2) = chi(2,0) = std::stod(value);
		}
		else if (param == "chi_yz")
		{
			// epsilon_scalar = epsilon0*std::stod(value);
			chi(1,2) = chi(2,1) = std::stod(value);
		}
		else if (param == "eps_vec")
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
							epsilon(k,k) = epsilon0*std::stod(coord);
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
		else if (param == "mu_vec")
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
							mu(k,k) = mu0*std::stod(coord);
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
		else if (param == "sigma_vec")
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
							sigma(k,k) = std::stod(coord);
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
		else if (param == "chi_vec")
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
							chi(k,k) = std::stod(coord);
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
			MyThrow(input_line,material_unknown_parameter);
	}

	//getters
	Eigen::Matrix3d Epsilon(void) { return epsilon; }
	Eigen::Matrix3d Mu(void) { return mu; }
	Eigen::Matrix3d Sigma(void) { return sigma; }
	Eigen::Matrix3d Chi(void) { return chi; }
	bool IsPML(void) { return pml; }
	
	private:
	double epsilon_scalar, sigma_scalar, mu_scalar, chi_scalar;
	bool pml;
	Eigen::Matrix3d epsilon, sigma, mu, chi;
};
