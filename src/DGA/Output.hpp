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
 
 //file Simulation.hpp
 #include "Utilities.hpp"

class Output
{
	public:
	Output()
	{
		// probepoint = Eigen::Vector3d({0,0,0});
		mode = "none";
		name = "simulation";
		output_period = 0;
		index = 0;
	}
	
	void Initialize(void)
	{
		index = 0;
	}
	
	const uint32_t Nprobes(void) const { return probepoint.size(); }
	
	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "mode")
		{
			if (std::find(outputmodes.begin(),outputmodes.end(),value) == outputmodes.end())
				MyThrow(input_line,sim_unknown_output);
			mode = value;
		}
		else if (param == "radiator")
			radiating_vol_bnd.push_back(std::stod(value));
		else if (param == "probe")
		{
			Eigen::Vector3d new_probepoint;
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
							new_probepoint[k]= std::stod(coord);
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
			
			probepoint.push_back(new_probepoint);
		}
		else if (param == "period")
		{
			output_period = std::stod(value);
		}
		else if (param == "name")
			name = value;
		else
			MyThrow(input_line,out_unknown_parameter);
	}
	
	bool AllowPrint(double t)
	{
		if (t >= index)
		{
			index+=output_period;
			return true;
		}
		else
			return false;
	}
	
	const double& Period(void) const { return output_period; }
	const std::string& Name(void) const { return name; }
	const OutputMode& Mode(void) const { return mode; }
	const std::vector<uint32_t> GetRadiators(void) const {return radiating_vol_bnd; };
	const Eigen::Vector3d& Probe(uint32_t i) const
	{ 
		assert(i<probepoint.size());
		return probepoint[i];
	}
	
	private:
	std::vector<uint32_t> radiating_vol_bnd;
	std::string name;
	double index;
	double output_period;
	std::vector<Eigen::Vector3d> probepoint;
	OutputMode mode;
};
