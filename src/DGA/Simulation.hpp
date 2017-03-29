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

class Simulation
{
	public:
	Simulation()
	{
		d=0;
		sources=std::vector<uint32_t>({1});
		mesh_label=1;
		output = 0;
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
			output = std::stod(value);
		// else if (param == "probe")
		// {
			// auto  i = value.begin();
			// if (*i != '{')
				// MyThrow(input_line,coordinates_syntax);
			// else
			// {
				// uint8_t k=0;
				// i++;
				// while (*i != ',' && *i != '}' && i != value.end())
				// {
					// std::string coord;
					// while (*i != ',' && *i != '}' && i != value.end())
					// {
						// coord.push_back(*i);
						// i++;
					// }
					// if (i == value.end())
						// MyThrow(input_line,unbalanced_bracket);
					// else 
					// {
						// if (k < 3)
						// {
							// probepoint[k]= std::stod(coord);
							// k++;
						// }
						// else
							// MyThrow(input_line,too_many_coords);
						
						// if (*i == ',')
							// i++;
					// }
					
				// }
				
				// if (i == value.end())
					// MyThrow(input_line,too_few_coords);
			// }
		// }
		else
			MyThrow(input_line,sim_unknown_parameter);
	}
	
	const Duration& Time(void) const { return d; } 
	const uint32_t& MeshLabel(void) const { return mesh_label; }
	const uint32_t& Output(void) const { return output; }
	const std::vector<uint32_t>& Src(void) { return sources; }
	// const Eigen::Vector3d& Probe(void) const { return probepoint; }
	
	private:
	Duration d;
	// Eigen::Vector3d probepoint;
	std::vector<uint32_t> sources; //can combine multiple sources
	uint32_t mesh_label, output;
};
