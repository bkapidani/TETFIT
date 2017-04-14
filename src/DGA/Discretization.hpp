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
#include "Utilities.hpp" //contains also all includes to c++ std libraries, EIGEN and SILOs
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "Simulation.hpp"
#include "Output.hpp"
#include "bessel.h"

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
								else if (thing_being_defined == "output")
									Outputs[definition_label]=Output();
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
						else if (thing_being_defined == "output")
							Outputs[definition_label].SetParam(input_line,tok,val);
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
		
		Run();
	}

	void Run(void)
	{
		timecounter t_sim;
		t_sim.tic();
		
		for (auto sims : Simulations)
			RunSimulation(sims.second, sims.first);
		t_sim.toc();
		
		std::cout << std::endl;
		std::cout << "----------- All simulations ran. Elapsed time: " << t_sim << " seconds ----------" << std::endl;
	}
	
	double analytic_value(SpaceTimePoint p, double sigma, double eps, double mu, double freq)
	{
		auto x = p[0]; auto y = p[1]; auto z = p[2]; auto t = p[3];
		double c = 1/sqrt(eps*mu);
		double ax=5e-2;
		double az=10e-2;
		double ay=2.5e-2; //momentarily useless
		double ksi = sigma/eps/2;
		double alph1 = pow(3.141592*c/ax,2);
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
		
		int32_t lim1 = floor( (t * c - z) / (2 * az) );
		int32_t lim2 = floor( (t * c + z) / (2 * az) - 1);
		
		// std::cout << "At time t = " << t << " limits are : " << lim1 << "\t" << lim2 << std::endl;
		
		double a1, a2;
		int32_t j=0;
		// i=j=0;
		a1=a2=0;
		
		
		while (true)
		{
			double k = (z + 2 * az * j) / c;
			if (t>k)
			{
				// std::cout << "First set of waves: integrate from " << k << " to " << t << " seconds" << std::endl;
				// std::cout << "Alpha =" << alpha << std::endl;
				a1 += sin(3.141592*x/ax)*inverse_laplace_transform(t,k, alpha, ksi, c, freq,flag);
			}
			else
				break;
			j++;
		}
		
		j=1;
		while (true)
		{
			double k = (2*az*j - z) / c;
			if (t>k)
			{
				// std::cout << "Second set of waves: integrate from " << k << " to " << t << " seconds" << std::endl;
				// std::cout << "Alpha =" << alpha << std::endl;
				a2 -= sin(3.141592*x/ax)*inverse_laplace_transform(t,k, alpha, ksi, c, freq,flag);
			}
			else
				break;
			j++;
		}
		
		return a1+a2;
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
				std::vector<std::vector<uint32_t>>().swap(edge_src);
				std::vector<std::vector<uint32_t>>().swap(face_src);
				std::vector<double>().swap(CellVolumes);
				// std::vector<double>().swap(probe_numeric_times);
				// std::vector<double>().swap(probe_numeric_yvalues);
				// std::vector<double>().swap(probe_numeric_zvalues); 
				// std::vector<double>().swap(probe_numeric_xvalues);
				std::vector<Eigen::Vector3d>().swap(edge_bars);
				std::vector<Eigen::Vector3d>().swap(face_bars);
			}
			else
			{
				//Still need to merge variables
			}
			m.Switch();
		}
	}
	
	void RunSimulation(const Simulation& s, const uint32_t& sim_label)
	{
		std::cout  << std::endl << std::endl; 
		std::cout <<"------------------------ Running FDTD simulation ------------------------" << std::endl << std::endl;
		
		current_simulation = sim_label;
		// auto sim_sources = s.Src();
		auto m = Meshes[s.MeshLabel()];
		auto o = &Outputs[s.Output()];
		auto mod_out = (*o).Mode();
		(*o).Initialize();
		Duration simulation_time = s.Time();
		
		meshlock.lock(); //lock access to the meshes map
		
		FlushMesh();
		ReadMesh(m);
		loaded_mesh_label = s.MeshLabel();
		t_step = estimate_time_step_bound();
		ConstructMaterialMatrices();
		
		// Eigen::VectorXd rhs_PM = Eigen::VectorXd::Random(edges_size());
		// double new_norm = std::pow(double(2)/t_step,double(2))/(rhs_PM.norm());
		// rhs_PM = new_norm*rhs_PM;
		timecounter t_power;
		// t_power.tic();
		// double t_step_power = estimate_time_step_bound_algebraic(rhs_PM);
		// t_power.toc();
		// std::cout << "Power method took " << t_power << "s." << std::endl;
		t_power.tic();
		double t_step_spectra = estimate_ts_w_spectra();
		t_power.toc();
		std::cout << "Spectral method took " << t_power << "s." << std::endl;
		
		auto t_step_geom = t_step;
		t_step = 0.98*t_step_spectra;
		ConstructMaterialMatrices();
		
		if ( mod_out == "probepoint")
		{
			// std::cout << "ci sono" << std::endl;
			probe_numeric_times.clear();
			probe_elem.clear();
			probe_numeric_xvalues.resize(0);
			probe_numeric_yvalues.resize(0);
			probe_numeric_zvalues.resize(0);
			analytic_value_vector.resize(0);
			probe_points.resize(0);
			
			std::vector<double> dummy_probe;
			
			for (uint32_t p=0; p<(*o).Nprobes(); p++)
			{
				auto pv = (*o).Probe(p);
				probe_points.push_back(pv);
				probe_elem.push_back(FindProbe(pv));
				probe_numeric_xvalues.push_back(dummy_probe);
				probe_numeric_yvalues.push_back(dummy_probe);
				probe_numeric_zvalues.push_back(dummy_probe);
				analytic_value_vector.push_back(dummy_probe);
				
			}
			
			// std::cout << "ci sono" << std::endl;
		}
		
		// for (auto bb : BCs)
			// std::cout << "Boundary condition: " << bb.first << " of type " << bb.second.Type() << std::endl;
		// for (auto bb : Sources)
			// std::cout << "Source: " << bb.second.Surface() << " of type " << bb.second.Type() << std::endl;
		
		double step_time_average=0;
		const uint32_t N_of_steps=simulation_time/t_step;
		uint64_t i;
		
		std::cout << std::endl << "Simulation parameters:" 		<< std::endl;
		std::cout << std::setw(20) << "Mesh: "             		<< std::setw(20) << m.FileName()              				  << std::endl;
		std::cout << std::setw(20) << "Minimum diameter: " 		<< std::setw(20)  << min_h 						<< "   m" << std::endl;
		std::cout << std::setw(20) << "Simulation time: "  		<< std::setw(20) << simulation_time              << " sec" << std::endl;
		std::cout << std::setw(20) << "Time step geom: "   		<< std::setw(20) << t_step_geom                      << " sec" << std::endl;
		// std::cout << std::setw(20) << "Time step power: "  		<< std::setw(20) << t_step_power                       << " sec" << std::endl;
		std::cout << std::setw(20) << "Time step spectral: "  	<< std::setw(20) << t_step                       << " sec" << std::endl;
		std::cout << std::setw(20) << "Unknowns: "         		<< std::setw(20) << U_frac_size+F_frac_size      << std::endl  << std::endl;
		
		
		// T time_function;
		timecounter step_cost;
		std::vector<double> numeric_values,numeric_times;
		U = I = Eigen::VectorXd::Zero(edges_size());
		F = Eigen::VectorXd::Zero(surfaces_size());
		Eigen::VectorXd curl_u(surfaces_size()), curl_f(edges_size());
		Eigen::VectorXd U_frac = Eigen::VectorXd::Zero(U_frac_size);
		Eigen::VectorXd F_frac = Eigen::VectorXd::Zero(F_frac_size);
		
		auto Old_F_frac = F_frac;
		auto Old_U_frac = U_frac;
		
		auto start_of_u = U_frac.data(); //pointer to the start of the big electric vector
		auto start_of_f = F_frac.data(); //pointer to the start of the big magnetic vector
		
		Eigen::Map<Eigen::VectorXd> U_a(start_of_u,H_size), U_b(start_of_u+H_size,P_size), U_c(start_of_u+H_size+P_size,Q_size);
		Eigen::Map<Eigen::VectorXd> U_d(start_of_u+H_size+P_size+Q_size,B_size);
		Eigen::Map<Eigen::VectorXd> F_a(start_of_f,N_size), F_b(start_of_f+N_size,R_size), F_c(start_of_f+N_size+R_size,S_size);
		
		double current_time;
		// auto pippo = M.transpos()*C.transpose();
		// std::ofstream file_pippo("Stau_per_C.dat");
		// filo_pippo << pippo;
		// file_pippo.close();
		// std::cout << M.norm() << " " << edges.size() << std::endl;
		// Actual simulation!
		for (i=1; i*t_step <= simulation_time; i++)
		{
			step_cost.tic();
			current_time = double(i)*t_step;
			
			// std::cout << U_a.lpNorm<Eigen::Infinity>() << '\t' << U_b.lpNorm<Eigen::Infinity>() << '\t'
			          // << U_c.lpNorm<Eigen::Infinity>() << '\t' << U_d.lpNorm<Eigen::Infinity>() << std::endl;
			
			for (uint32_t ee = 0; ee < edges_size(); ee++)
			{
				if (edge_bcs[ee] != 0  && BCs[edge_bcs[ee]].Type() != "none")
				{
					U_d[boundary_index[ee]] = ComputeEdgeBC(ee,current_time); //only boundary edges can have boundary conditions
					
					// for (auto ii : associated_h_edges[ee])
						// U_a[ii] = ubc;
					// for (auto ii : associated_p_edges[ee])
						// U_b[ii] = ubc;
					// for (auto ii : associated_frac_edges[ee])
						// U_c[ii] = 0.5*ubc;
				}
				else if (edge_src[ee].size()>0)
				{
					auto isrc = ComputeCurrentSource(ee,current_time - 0.5*t_step);
					auto usrc = ComputeEfieldSource(ee, current_time - t_step);
					I[ee] = isrc;
					
					for (auto ii : associated_h_edges[ee])
					{
						U_a[ii] = usrc;
						// std::cout << "debuggy buggy" << std::endl;
					}
					for (auto ii : associated_p_edges[ee])
					{
						U_b[ii] = usrc;
						// std::cout << "debuggy buggy" << std::endl;
					}
					for (auto ii : associated_frac_edges[ee])
					{
						U_c[ii] = 0.5*usrc;
						// std::cout << "debuggy buggy" << std::endl;
					}
					for (auto ii : associated_bnd_edges[ee])
					{
						U_d[ii] = usrc;
						// std::cout << "debuggy not buggy" << std::endl;
					}
				}
			}
			
			// std::cout << U_a.lpNorm<Eigen::Infinity>() << '\t' << U_b.lpNorm<Eigen::Infinity>() << '\t'
			          // << U_c.lpNorm<Eigen::Infinity>() << '\t' << U_d.lpNorm<Eigen::Infinity>() << std::endl;		  
			// std::cout << std::endl;
			
			// std::cout << U_frac.lpNorm<Eigen::Infinity>() << "\t";
			U = M*U_frac;
			// assert(U_frac.lpNorm<Eigen::Infinity>() == U.lpNorm<Eigen::Infinity>());
			// std::cout << U.lpNorm<Eigen::Infinity>() << std::endl;
			
			Old_F_frac = F_frac;
			
			// Magnetic Part:
			curl_u     = C*U;
			
			F_a       -= t_step*N*curl_u;
			F_b        = R*Old_F_frac - t_step*Tr*curl_u;
			F_c        = S*F_c - t_step*Ts*curl_u;
			F          = T*F_frac;
			
			Old_U_frac = U_frac;
			
			// Electric Part:			
			curl_f     = C.transpose()*F-I;
			
			U_a       += t_step*H*curl_f;
			U_b        = P_p.cwiseProduct(U_b) + t_step*Mp*curl_f;
			U_c        = Q*U_c + t_step*Mq*curl_f;
			
			step_cost.toc();
			step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();

			// Debug
			// std::cout << "Time: "      << std::setw(20) << current_time << '\t'; 
			// std::cout << "Maximum F: " << std::setw(20) << F.lpNorm<Eigen::Infinity>() << '\t'; 
			// std::cout << "Maximum U: " << std::setw(20) << U.lpNorm<Eigen::Infinity>() << std::endl;
			
			if ((*o).AllowPrint(current_time))
				ExportFields(mod_out, current_time);

			if (i % 140 == 0)
				std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(9) << step_time_average << "s, " 
						  << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
		}
		
		meshlock.unlock(); //unlock the access to the meshes map

		// /* Output stats and fields*/
		if (mod_out == "probepoint")
		{
			std::ofstream os, os_a;
			
			for (uint32_t p=0; p<probe_elem.size(); p++)
			{
				std::stringstream ss, sa;
				ss << (*o).Name() << std::setw(5) << std::setfill('0') << sim_label << "_probe" << p+1 << ".dat";
				sa << (*o).Name() << std::setw(5) << std::setfill('0') << sim_label << "_probe" << p+1 << "_analytic.dat";
				os.open(ss.str().c_str());
				os_a.open(sa.str().c_str());
				
				for (size_t k=0; k < probe_numeric_times.size(); k++)
				{
				  os << std::setw(15) << probe_numeric_times[k];
				  os << std::setw(15) << probe_numeric_xvalues[p][k];
				  os << std::setw(15) << probe_numeric_yvalues[p][k];
				  os << std::setw(15) << probe_numeric_zvalues[p][k];
				  os << std::endl;
				  
				  os_a << std::setw(15) << probe_numeric_times[k] << " " << std::setw(15) << analytic_value_vector[p][k] <<  std::endl;
				}
				
				os.close();
				os_a.close();
			}
		}
		// else if (mod_out == "silo")
		// {

		// }
		
		std::cout << std::setw(20) << "Time step average cost is "  << step_time_average/(double(i)) 
		          << " seconds (" << i << " time steps!)" << std::endl;
		std::cout << std::setw(20) << "Total running time is "      << step_time_average << " seconds" << std::endl;
	}
	
    bool ExportMesh(std::string meshname)
    {
        if (!_siloDb)
        {
            std::cout << "Silo database not opened" << std::endl;
            return false;
        }
        
		
        if (meshname.size() == 0)
        {
            std::cout << "Mesh without name, cannot export!" << std::endl;
            return false;
        }
        
        using namespace std::chrono;
        high_resolution_clock::time_point start, stop;
        duration<double> time_span;
        
        start = high_resolution_clock::now();
        
        /* Step 1: Make node arrays.
         * 1 -> x1,y1,z1
         * 2 -> x2,y2,z2
         * ...
         */
        std::vector<double> coords_x, coords_y, coords_z;
        
        coords_x.reserve( pts.size() );
        coords_y.reserve( pts.size() );
        coords_z.reserve( pts.size() );
        
        for (auto itor = pts.begin(); itor != pts.end(); itor++)
        {
            coords_x.push_back( (*itor)(0) );
            coords_y.push_back( (*itor)(1) );
            coords_z.push_back( (*itor)(2) );
        }
        
        /* Step 2: Make nodelist: |-t1-|-t2-|....|-tn-| where -tn- is the
         * quadruple of nodes composing the n-th tetrahedron.
         */
        
        std::vector<int> nodelist;
        nodelist.reserve( volumes_size() );
        
        for (auto itor = 0; itor < volumes_size(); itor++)
        {
            auto vol = volumes[itor];
            auto ptids = std::vector<uint32_t>({std::get<0>(vol),std::get<1>(vol),std::get<2>(vol),std::get<3>(vol)});
            
            if ( CellVolumes[itor] > 0 )
            {
                nodelist.push_back( static_cast<int>(ptids.at(0)) + 1 );
                nodelist.push_back( static_cast<int>(ptids.at(1)) + 1 );
                nodelist.push_back( static_cast<int>(ptids.at(2)) + 1 );
                nodelist.push_back( static_cast<int>(ptids.at(3)) + 1 );
            }
            else
            {
                nodelist.push_back( static_cast<int>(ptids.at(0)) + 1 );
                nodelist.push_back( static_cast<int>(ptids.at(2)) + 1 );
                nodelist.push_back( static_cast<int>(ptids.at(1)) + 1 );
                nodelist.push_back( static_cast<int>(ptids.at(3)) + 1 );
            }
        }
        
        /* Step 3: Put mesh data in the database. */
        int shapesize = 4;
        int shapecounts = static_cast<int>( volumes_size() );
        int nnodes = static_cast<int>( pts.size() );
        int nzones = static_cast<int>( volumes_size() );
        int ndims = 3;
        
        DBPutZonelist(_siloDb, "zonelist", nzones, ndims, nodelist.data(),
                      static_cast<int>( nodelist.size() ), 1,
                      &shapesize, &shapecounts, 1);
        
		double *coords[] = { coords_x.data(), coords_y.data(), coords_z.data() };
		
		DBPutUcdmesh(_siloDb, meshname.c_str(), ndims, NULL, coords,
					 nnodes, nzones, "zonelist", NULL, DB_DOUBLE, NULL);
        
        stop = high_resolution_clock::now();
        time_span = duration_cast<duration<double>>(stop - start);
        
        std::cout << "SILO: Mesh export done in " << time_span.count() << " seconds.";
        std::cout << std::endl;
        
        return true;
    }

	void ExportFields(const std::string s, double t)
	{
		// Eigen::Vector3d num_val;
		if (s == "probepoint")
		{
			for (uint32_t p=0; p<probe_elem.size(); p++)
			{
				auto num_ele = GetElectricField(probe_elem[p]);
				probe_numeric_xvalues[p].push_back(num_ele[0]);
				probe_numeric_yvalues[p].push_back(num_ele[1]);
				probe_numeric_zvalues[p].push_back(num_ele[2]);
				SpaceTimePoint stp = SpaceTimePoint({probe_points[p][0],probe_points[p][1],probe_points[p][2],t});
				double anal_value = analytic_value(stp,Materials[vol_material[probe_elem[p]]].Sigma(),
				                                       Materials[vol_material[probe_elem[p]]].Epsilon(),
													   Materials[vol_material[probe_elem[p]]].Mu(),5e9); //BIG HACK!
				analytic_value_vector[p].push_back(anal_value);
			}
			
			probe_numeric_times.push_back(t);
		}
		else if (s == "silo")
		{
			timecounter t_export;
			t_export.tic();
			
			auto meshname = Meshes[loaded_mesh_label].Name();
			auto op = Outputs[Simulations[current_simulation].Output()];
			
			std::stringstream ss;
			ss << op.Name() 
			   << std::setw(5) << std::setfill('0') << current_simulation << "_" 
			   << std::setw(5) << std::setfill('0') << uint32_t(t/t_step) << ".silo";
			
			auto filename = ss.str();
			_siloDb = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
			
			if (!ExportMesh(meshname))
				std::cout << "Problems with mesh export!" << std::endl;
			
			std::vector<double> Ex_vals, Ey_vals, Ez_vals, Hx_vals, Hy_vals, Hz_vals;
			Ex_vals.reserve( volumes_size() );
			Ey_vals.reserve( volumes_size() );
			Ez_vals.reserve( volumes_size() );
			// Hx_vals.reserve( volumes_size() );
			// Hy_vals.reserve( volumes_size() );
			// Hz_vals.reserve( volumes_size() );
			
			for (auto itor = 0; itor < volumes_size(); itor++)
			{
				auto Efield = GetElectricField(itor);
				// auto Hfield = GetMagneticField(itor);
				Ex_vals.push_back(Efield(0));
				Ey_vals.push_back(Efield(1));
				Ez_vals.push_back(Efield(2));
			}
			
			std::vector<std::string> varnames({"Ex","Ey","Ez","Hx","Hy","Hz"});
			
			DBPutUcdvar1(_siloDb, varnames[0].c_str(), meshname.c_str(), Ex_vals.data(), static_cast<int>(volumes_size()), NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
			DBPutUcdvar1(_siloDb, varnames[1].c_str(), meshname.c_str(), Ey_vals.data(), static_cast<int>(volumes_size()), NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
			DBPutUcdvar1(_siloDb, varnames[2].c_str(), meshname.c_str(), Ez_vals.data(), static_cast<int>(volumes_size()), NULL, 0, DB_DOUBLE, DB_ZONECENT, NULL);
			
			DBClose(_siloDb);
			
			t_export.toc();
			std::cout << "SILO: Output to file done in " << std::setw(7) << t_export << std::setw(8) << " seconds" << std::endl;
		}
	}
	
	uint32_t FindProbe(const Eigen::Vector3d& probe_vec)
	{
		double dist = (vol_barycenter(0)-probe_vec).norm();
		uint32_t ret = 0;
		
		for (uint32_t v=0; v<volumes_size(); v++)
		{
			double new_dist = (vol_barycenter(v)-probe_vec).norm();
			
			if (new_dist < dist)
			{
				dist = new_dist;
				ret = v;
			}
		}

		// std::ofstream debug_faces("debug_faces.txt", std::ofstream::out | std::ofstream::app);
		
		// for (auto ff : vtf_list[ret])
			// debug_faces << print_face(4,abs(ff),true,0,255,255);
		
		// debug_faces.close();
		
		return ret;
	}
	
	Eigen::Vector3d GetElectricField(uint32_t vol)
	{	
		std::array<uint32_t,6>       edgs = primal_edge_labels(vol);
		std::array<Eigen::Vector3d,4> pav = primal_area_vectors(vol);
		double v = 3*CellVolumes[vol];
		Eigen::Vector3d ret(0,0,0);
		
		std::vector<double> eq;
		for (auto e : edgs)
			eq.push_back(U[e]);
		
		ret += ( pav[1]*eq[0] - pav[2]*eq[2] + pav[3]*eq[3])/v;
		ret += ( pav[0]*eq[0] - pav[2]*eq[1] + pav[3]*eq[4])/v;
		ret += (-pav[1]*eq[1] + pav[0]*eq[2] + pav[3]*eq[5])/v;
		ret += ( pav[0]*eq[3] - pav[1]*eq[4] + pav[2]*eq[5])/v;
		
		return ret/4;
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
		uint32_t counter = 0;
		
		for (auto src : edge_src[e])
		{
			counter++;
			if (Sources[src].Type() == "e")
				vector_val += Sources[src].Compute(p);
		}
		
		// std::cout << "Visited loop " << counter << " times!" << std::endl;
		
		// std::cout << "-------------" << std::endl << vector_val[1] << std::endl;
		
		return vector_val.dot(pts[abs(etn_list[e][1])]-pts[abs(etn_list[e][0])]);
	}

	double ComputeCurrentSource(uint32_t e, double t)
	{
		auto eb = edge_barycenter(e);
		SpaceTimePoint p;
		p[0] = eb[0];
		p[1] = eb[1];
		p[2] = eb[2];
		p[3] = t;
		
		Eigen::Vector3d vector_val(0,0,0);
		uint32_t counter = 0;
		
		for (auto src : edge_src[e])
		{
			counter++;
			if (Sources[src].Type() == "j")
				vector_val += Sources[src].Compute(p);
		}
		
		// std::cout << "Visited loop " << counter << " times!" << std::endl;
		
		// std::cout << "-------------" << std::endl << vector_val[1] << std::endl;
		
		return vector_val.dot(pts[abs(etn_list[e][1])]-pts[abs(etn_list[e][0])]); //wrong, integrate on dual faces
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
		
		auto ee = arr[0].first;
		uint32_t       p0(std::get<0>(ee));
		uint32_t       p1(std::get<1>(ee));
		sgnint32_t<int32_t> n1(p0,-1);
		sgnint32_t<int32_t> n2(p1, 1);
		
		std::vector<sgnint32_t<int32_t>> dummy(2);		
		etn_list.push_back(dummy);
		etn_list[edges.size()][0] = n1;
		etn_list[edges.size()][1] = n2;
		sgnint32_t<int32_t> e1(edges.size(),-1);
		sgnint32_t<int32_t> e2(edges.size(), 1);
		dual_star_offsets.push_back(std::vector<uint32_t>({uint32_t(nte_list[p0].size()),uint32_t(nte_list[p1].size())}));
		nte_list[p0].push_back(e1);
		nte_list[p1].push_back(e2);
		edges.push_back(ee);
		
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

				auto ee = arr[left].first;
				uint32_t       p0(std::get<0>(ee));
				uint32_t       p1(std::get<1>(ee));
				sgnint32_t<int32_t> n1(p0,-1);
				sgnint32_t<int32_t> n2(p1, 1);
				
				std::vector<sgnint32_t<int32_t>> dummy(2);		
				etn_list.push_back(dummy);
				etn_list[edges.size()][0] = n1;
				etn_list[edges.size()][1] = n2;
				sgnint32_t<int32_t> e1(edges.size(),-1);
				sgnint32_t<int32_t> e2(edges.size(), 1);
				
				dual_star_offsets.push_back(std::vector<uint32_t>({uint32_t(nte_list[p0].size()),uint32_t(nte_list[p1].size())}));
				nte_list[p0].push_back(e1);
				nte_list[p1].push_back(e2);
				
				edges.push_back(ee);
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
		{
			ebs[i] =  edge_barycenter(edgs[i]);
			// std::cout << "(" << ebs[i][0] << ", " << ebs[i][1] << ", " << ebs[i][2] << ")" << std::endl;
		}
		
		// std::cout << "-------" << std::endl;
		
		std::array<Eigen::Vector3d, 4> fbs;
		for (uint32_t i = 0; i < 4; i++)
		{
			fbs[i] = face_barycenter(fcs[i]);
			// std::cout << "(" << fbs[i][0] << ", " << fbs[i][1] << ", " << fbs[i][2] << ")" << std::endl;
		}

		
		Eigen::Vector3d vb = vol_barycenter(vol);
		
		/* Area of quadrilateral ABCD is 0.5*|AC x BD|. Orientation of the dual
		 * area vector must be the same of the primal edge vector i.e. it must
		 * satisfy 'dot(PEV, DAV) >= 0' */
		
		double sign = CellVolumes[vol]<0 ? -0.5 : 0.5;
		ret[0] = ((vb-ebs[0]).cross(fbs[2]-fbs[3]))*sign;
		ret[1] = ((vb-ebs[1]).cross(fbs[0]-fbs[3]))*sign;
		ret[2] = ((vb-ebs[2]).cross(fbs[3]-fbs[1]))*sign;
		ret[3] = ((vb-ebs[3]).cross(fbs[1]-fbs[2]))*sign;
		ret[4] = ((vb-ebs[4]).cross(fbs[2]-fbs[0]))*sign;
		ret[5] = ((vb-ebs[5]).cross(fbs[0]-fbs[1]))*sign;
		
		return ret;
	}

	std::array<Eigen::Vector3d, 6> primal_edge_vectors(const uint32_t& vol)
	{
		std::vector<uint32_t> edgs;
		std::array<Eigen::Vector3d, 6> ret;
		
		for (auto ff : vtf_list[vol])
		{
			
			for (auto ee : fte_list[abs(ff)])
			{
				edgs.push_back(abs(ee));
			}
		}
		sort_unique(edgs);
		std::swap(edgs[1],edgs[2]);
		std::swap(edgs[1],edgs[3]);
		
		for (uint32_t i=0; i<6; i++)
			ret[i]=pts[abs(etn_list[edgs[i]][1])]-pts[abs(etn_list[edgs[i]][0])];
		return ret;
	}
	
	std::array<uint32_t, 6> primal_edge_labels(const uint32_t& vol)
	{
		std::vector<uint32_t> edgs;
		std::array<uint32_t, 6> ret;
		
		for (auto ff : vtf_list[vol])
		{
			
			for (auto ee : fte_list[abs(ff)])
			{
				edgs.push_back(abs(ee));
			}
		}
		sort_unique(edgs);
		std::swap(edgs[1],edgs[2]);
		std::swap(edgs[1],edgs[3]);
		
		for (uint32_t i=0; i<6; i++)
			ret[i]=edgs[i];
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
		if (msh.FileName().size() == 0)
		{
			std::cout << "Invalid mesh file name" << std::endl;
			return false;
		}
		
		uint32_t	lines, linecount;
		
		mapped_file mf(msh.FileName());
		
		// std::cout << " * * * Reading NETGEN format mesh * * * " << std::endl;
		
		tctot.tic();
		
		/************************ Read points ************************/
		linecount = 0;
		
		const char *data = mf.mem();
		char *endptr;
		
		lines = strtot<uint32_t>(data, &endptr);
		std::vector<Eigen::Vector3d> new_pts;
		pts = std::move(new_pts);
		pts.reserve(lines);
		std::vector<std::vector<uint32_t>> ass_vols;
		this->associated_volumes = std::move(ass_vols);
		
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
		
		std::vector<volume_type> vlms;
		volumes = std::move(vlms);
		volumes.reserve(lines);
		std::vector<std::pair<double,double>> dual_vol_parameters(pts.size(),std::make_pair<double,double>(-1,-1));
		std::bitset<3> bob;
		bob[0]=0; bob[1]=0; bob[2]=0;
		std::vector<std::bitset<3>> classify_nodes(pts.size(),bob);
		std::vector<bool> bnd_nodes(pts.size(),false);
		
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
		
		std::vector<std::pair<double,double>> primal_vol_parameters(linecount,std::make_pair<double,double>(-1,-1));
		
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
			
			auto mat_label = std::get<1>(tet);
			vol_material.push_back(mat_label);
			
			auto pp = std::vector<uint32_t>({p0,p1,p2,p3});
			
			for (uint8_t ip=0; ip<4; ip++)
			{
				if (dual_vol_parameters[pp[ip]].first!=-1)
				{
					if (Materials[mat_label].Sigma()   != dual_vol_parameters[pp[ip]].second)
						classify_nodes[pp[ip]][2] = true;
					if (Materials[mat_label].Epsilon() != dual_vol_parameters[pp[ip]].first)
						classify_nodes[pp[ip]][1] = true;
					if (dual_vol_parameters[pp[ip]].second != 0)
					{
						classify_nodes[pp[ip]][0] = true;
					}
				}
				else
				{
					dual_vol_parameters[pp[ip]].first  = Materials[mat_label].Epsilon();
					dual_vol_parameters[pp[ip]].second = Materials[mat_label].Sigma();
					if (dual_vol_parameters[pp[ip]].second != 0)
						classify_nodes[pp[ip]][0] = true;
				}
				
			}
		}
		
		//std::cout << volumes_size() << std::endl;
		
		std::vector<tm_tuple>().swap(temp_tet);
		std::vector<uint32_t> labels(4*lines);
		std::vector<surface_type> srfcs;
		surfaces = std::move(srfcs);
		surfaces.reserve(4*lines);
		
		unique(temp_tri0, labels); //this also fills the surfaces vector
		std::vector<label_surface_type>().swap(temp_tri0);
		
		std::vector<cluster_list> dummy(surfaces_size()), dummyv;
		ftv_list = std::move(dummy);
		vtf_list = std::move(dummyv);
		vtf_list.reserve(volumes.size());
		
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
		std::vector<cluster_list> dummye;
		fte_list = std::move(dummye);
		fte_list.reserve(lines);
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
		
		std::vector<edge_type> dgs;
		edges = std::move(dgs);
		edges.reserve(3*lines);
		
		std::vector<cluster_list> dummyn2(pts.size());
		this->nte_list=std::move(dummyn2);
		std::vector<cluster_list> dummyn;
		etn_list = std::move(dummyn);
		
		std::vector<uint32_t> e_labels(3*lines);
		unique(temp_edge0, e_labels); //fills edges and etn, nte adjacency lists

		// etn_list.reserve(edges_size());
		
		std::vector<label_edge_type>().swap(temp_edge0);
		std::vector<cluster_list> dummye2(edges_size());
		etf_list = std::move(dummye2);
		
		std::vector<double_triplet> tripletList;
		tripletList.reserve(3*surfaces_size());
		Eigen::SparseMatrix<double> C(surfaces_size(),edges_size());
		n_index.resize(surfaces_size(),0);
		r_index = n_index;
		// s_index.resize(surfaces_size());
		N_size = R_size = S_size = 0;
		std::vector<uint32_t> dual_is_fractured(nodes_size(),0), primal_is_fractured(volumes_size(),0);
		std::vector<bool> bnd_edges(edges.size(),false);
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
			// std::cout << "k is " << k << " and surfaces.size() is " << surfaces.size() << std::endl;
			
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
					auto srf = surfaces[k];
					bnd_edges[abs(e1)]=true;
					bnd_edges[abs(e2)]=true;
					bnd_edges[abs(e3)]=true;
					bnd_nodes[std::get<0>(srf)]=true;
					bnd_nodes[std::get<1>(srf)]=true;
					bnd_nodes[std::get<2>(srf)]=true;
					
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
				
				if (!primal_is_fractured[vol1])
				{
					primal_is_fractured[vol1]=S_size+1;
					S_size+=4;
				}

				if (vols.size()==2)
				{
					auto vol2= abs(vols[1]);
					
					if (!primal_is_fractured[vol2])
					{
						primal_is_fractured[vol2]=S_size+1;
						S_size+=4;
					}
				
				}
			}
		}
		
		std::vector<uint32_t>().swap(e_labels);
		h_index.resize(edges.size(),0);
		p_index= boundary_index = h_index;
		// q_index.resize(edges.size());
		
		C.setFromTriplets(tripletList.begin(),tripletList.end());
		this->C=std::move(C);
		std::vector<double_triplet>().swap(tripletList);
		std::vector<std::vector<uint32_t>> associated_frac_edges(edges.size()), associated_p_edges(edges.size()); 
		std::vector<std::vector<uint32_t>> associated_bnd_edges(edges.size()),  associated_h_edges(edges.size());
		H_size = Q_size = P_size = B_size = 0;
		
		// domain splitting
		for (uint32_t k = 0; k < edges.size(); k++)
		{
			auto t = edges[k];
			
			uint32_t       p0(std::get<0>(t));
			uint32_t       p1(std::get<1>(t));
			
			auto rec0 = classify_nodes[p0].to_ulong();
			auto rec1 = classify_nodes[p1].to_ulong();
			
			if (bnd_edges[k])
			{
				classify_edges.push_back(4);
				boundary_index[k] = B_size++;
				associated_bnd_edges[k].push_back(boundary_index[k]);
			}
			else if (rec0>2 || rec1>2 || (rec1+rec0==1) || (rec1+rec0==3) 
				 || (dual_vol_parameters[p0].second != dual_vol_parameters[p1].second) )
			{
				classify_edges.push_back(3);

				if (dual_is_fractured[p0]==0)
				{
					dual_is_fractured[p0]=Q_size+1;
					auto n_star1 = nte_list[p0].size();
					associated_frac_edges[k].push_back(dual_star_offsets[k][0]+Q_size);
					Q_size+=n_star1;
				}
				else
					associated_frac_edges[k].push_back(dual_star_offsets[k][0]+(dual_is_fractured[p0]-1));

				if (dual_is_fractured[p1]==0)
				{
					dual_is_fractured[p1]=Q_size+1;
					auto n_star2 = nte_list[p1].size();
					associated_frac_edges[k].push_back(dual_star_offsets[k][1]+Q_size);
					Q_size+=n_star2;
				}
				else
					associated_frac_edges[k].push_back(dual_star_offsets[k][1]+(dual_is_fractured[p1]-1));
			}
			else if ( rec0 == 2 || rec1 == 2 || rec0 == 0 || rec1 == 0)
			{
				classify_edges.push_back(1);
				h_index[k]=H_size++;
				associated_h_edges[k].push_back(h_index[k]);
			}
			/*else if (bnd_nodes[p0] || bnd_nodes[p1])
			{
				classify_edges.push_back(3);

				if (dual_is_fractured[p0]==0)
				{
					dual_is_fractured[p0]=Q_size+1;
					auto n_star1 = nte_list[p0].size();
					associated_frac_edges[k].push_back(dual_star_offsets[k][0]+Q_size);
					Q_size+=n_star1;
				}
				else
					associated_frac_edges[k].push_back(dual_star_offsets[k][0]+(dual_is_fractured[p0]-1));

				if (dual_is_fractured[p1]==0)
				{
					dual_is_fractured[p1]=Q_size+1;
					auto n_star2 = nte_list[p1].size();
					associated_frac_edges[k].push_back(dual_star_offsets[k][1]+Q_size);
					Q_size+=n_star2;
				}
				else
					associated_frac_edges[k].push_back(dual_star_offsets[k][1]+(dual_is_fractured[p1]-1));
			}*/
			else /*if ( dual_vol_parameters[p0].second == dual_vol_parameters[p1].second )*/
			{
				classify_edges.push_back(2);
				p_index[k]=P_size++;
				associated_p_edges[k].push_back(p_index[k]);
			}

		}
		
		// std::cout << "# of boundary edges is: " << B_size << std::endl;
		
		// for ( auto nid : dual_is_fractured)
			// if (nid != 0)
				// std::cout << nid << " ";
		
		// std::cout << std::endl;
		
		this->associated_frac_edges = std::move(associated_frac_edges);
		this->associated_bnd_edges = std::move(associated_bnd_edges);
		this->associated_p_edges = std::move(associated_p_edges);
		this->associated_h_edges = std::move(associated_h_edges);
		this->dual_is_fractured   = std::move(dual_is_fractured);
		this->primal_is_fractured = std::move(primal_is_fractured);
		this->bnd_nodes = std::move(bnd_nodes);
		std::vector<uint32_t> dumb(edges_size(),0);
		// std::vector<bool> all_false(edges_size(),false);
		std::vector<bool> is_dirichlet(edges_size(),false);
		auto edge_sources_current_bid = dumb;
		// this->is_source = all_false;
		tc.toc();
		
		std::cout << "done - " << tc << " seconds" << std::endl;
		
		/************************ Read boundary surfaces ************************/
		linecount = 0;
		auto num_of_tets=lines;
		lines = strtot<uint32_t>(endptr, &endptr);
		
		face_bcs.resize(surfaces_size());
		edge_bcs.resize(edges_size());
		
		std::vector<std::vector<uint32_t>> dumb_edge(edges_size());
		std::vector<std::vector<uint32_t>> dumb_face(surfaces_size());
		
		face_src=dumb_face;
		edge_src=dumb_edge;
		edge_bids=dumb;
		// std::vector<bool> edge_done(edges_size(),false);
		uint32_t trinum=0;
		
		// std::ofstream debug_faces("debug_faces.txt");
		
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
				// debug_faces << print_face(1,face_label,true,0,255,0);
				for (auto ee : fte_list[face_label])
				{
					edge_bcs[abs(ee)] = bid;
					edge_src[abs(ee)].clear();
					is_dirichlet[abs(ee)] = true; //overrides sources
					// classify_edges[abs(ee)] = 4;
					// debug_faces << print_edge(3,abs(ee),true,255,0,255);
					// edge_done[abs(ee)]    = true; 
				}
			}
			else
			{	
				std::vector<uint32_t> edgs;
				for (auto ee : fte_list[face_label])
				{
					if (!is_dirichlet[abs(ee)])
					{
						edgs.push_back(abs(ee));
					}
				}
				
				for	(auto src_label : Simulations[current_simulation].Src())
				{
					auto src = Sources[src_label];
					if (src.Surface() == bid)
					{
						// debug_faces << print_face(2,face_label,true,255,0,0);
						if (src.Type() == "e" || src.Type() == "j")
						{
							for (auto ee : edgs)
							{
								if (std::find(edge_src[ee].begin(),edge_src[ee].end(),src_label) == edge_src[ee].end())
								{
									edge_src[ee].push_back(src_label);
									// classify_edges[ee]=4; //leave them alone like dirichlets
								}
							}
						}
					}
				}
			}

			
			linecount++;
		}
		
		tc.toc();
		
		// debug_faces.close();
		
		std::cout << "Reading triangles: " << linecount;
		std::cout << "/" << lines  << " - " << tc << " seconds"  << std::endl;
		
		tctot.toc();
		
		return true;
	}
	
	double estimate_time_step_bound() //theoretic bound
	{
		double ret=1;
		min_h=1e6;
		average_diameter=0;
		
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
				double diameter = 1.5*vol/(sqrt(ff(0)*ff(0) + ff(1)*ff(1) + ff(2)*ff(2)));
				double dt = diameter/c;
				
				if (dt<ret)
					ret=dt;
				if (diameter < min_h)
					min_h = diameter;
				
				average_diameter+=2*diameter;
			}			
		}
		
		// std::cout << "Minimum diameter: " << min_h << std::endl;
		min_h = 2*min_h;
		average_diameter /= (4*volumes_size());
		return ret;
	}
	
	Eigen::VectorXd power_method_iteration(Eigen::VectorXd b )
	{
		timecounter t_iter_pw;
		t_iter_pw.tic();
		
		
		
		// auto alfa = C*b;
		// auto beta = N*alfa;
		// auto gamma = C.transpose()*beta;
		// auto delta = Einv*gamma;
		
		auto delta = Einv*(C.transpose()*(N*(C*b)));
		
		t_iter_pw.toc();
		// std::cout << "iteration took " << t_iter_pw << "s." << std::endl;
 		return delta;
	}
	
	double estimate_ts_w_spectra(void) //use spectra
	{
		// We are going to calculate the eigenvalues of M
		// Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
		// Eigen::MatrixXd M = A + A.transpose();
		class MyMatProd
		{
			public:
			
			MyMatProd(uint32_t s, Eigen::SparseMatrix<double>& c, Eigen::SparseMatrix<double>& n, Eigen::SparseMatrix<double>& einv )
			{
				size = s;
				this->C = &c;
				this->N = &n;
				this->Einv = &einv;
			}
			
			int rows() { return size; }
			int cols() { return size; }
			
			void perform_op(double *x_in, double *y_out)
			{
				Eigen::Map<Eigen::VectorXd> X(x_in,size);
				Eigen::Map<Eigen::VectorXd> Y(y_out,size);
				
				Y = (*Einv)*((*C).transpose()*((*N)*((*C)*X)));
			}
			
			Eigen::SparseMatrix<double> *C, *N, *Einv;
			uint32_t size;
		};
		// Construct matrix operation object using the wrapper class MyMatProd
		MyMatProd op(this->edges_size(),C,N,Einv);

		// Construct eigen solver object, requesting the largest three eigenvalues
		Spectra::SymEigsSolver< double, Spectra::LARGEST_MAGN, MyMatProd > eigs(&op, 1, 10);

		// Initialize and compute
		eigs.init();
		int nconv = eigs.compute();

		// Retrieve results
		double lambda;
		if(eigs.info() == Spectra::SUCCESSFUL)
		{
			auto eig_vec = eigs.eigenvalues();
			lambda = eig_vec(0);
		}
		return double(2)/sqrt(lambda);
	}
	
	double estimate_time_step_bound_algebraic(Eigen::VectorXd b ) //power method
	{
		double lambda;
		double tol = 0.01;
		std::vector<double> values;
		uint32_t it = 0;
		

		Eigen::VectorXd bnew;

		while (true)
		{
			bnew = b/b.norm();
			b = Einv*(C.transpose()*(N*(C*bnew)));
			lambda = bnew.dot(b);
			
			values.push_back(double(2)/sqrt(lambda));
			
			if ((b-lambda*bnew).norm() <= tol*fabs(lambda))
				break;
		}
		
		std::ofstream dbg_pm("asymptote_power_method.dat");
		for (uint32_t i=0; i<values.size(); ++i)
			dbg_pm << std::setw(10) << i << " " << std::setw(20) << values[i] << std::endl;
		dbg_pm.close();
		double ts = double(2)/sqrt(lambda);
		return ts - tol*ts;
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
		std::vector<double_triplet> H_trip, Einv_trip, Mp_trip, Mq_trip, M_trip, P_trip, Q_trip;
		std::vector<double_triplet> N_trip, Tr_trip, Ts_trip, T_trip, R_trip, S_trip;  
		uint32_t jj,kk;
		
		auto local_mag_Id = Eigen::MatrixXd::Identity(4,4);
		for (uint32_t vv=0; vv < volumes.size(); vv++)
		{
			double mu_vol    = Materials[vol_material[vv]].Mu();
			double chi_vol	 = Materials[vol_material[vv]].Chi();
			uint32_t jj,kk;
			auto face_vecs = dual_area_vectors(vv);
			auto fids = vtf_list[vv];
			std::vector<uint32_t> abs_fids;			
			
			for (auto ff : fids)
				abs_fids.push_back(abs(ff));
			double coeff=mu_vol*36.0/(fabs(CellVolumes[vv]));
			
			local_M = local_Z = Eigen::Matrix4d::Zero();
			
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
				
			// auto local_N = (local_M + local_Z).inverse();
			Eigen::MatrixXd local_N = (local_M + local_Z).llt().solve(local_mag_Id);
			auto local_R = local_N*(local_M - local_Z);
			
			// std::cout << "Tetrahedron: [ " << std::get<0>(volumes[vv]) << " " << std::get<1>(volumes[vv]) << " " << std::get<2>(volumes[vv]) << " " << std::get<3>(volumes[vv]);
			// std::cout << " ]" << std::endl;
			// std::cout << std::endl << local_M << std::endl << std::endl;
			
			uint32_t offset;
			bool is_frac=false;
			if (primal_is_fractured[vv]>0)
			{
				is_frac=true;
				offset = primal_is_fractured[vv]-1;
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
							Ts_trip.push_back(double_triplet(offset+jj,abs_fids[kk],local_N(jj,kk)));
					}
				}
				
				jj++;
				kk=0;
			}
		}
		
		Eigen::VectorXd	P_p(P_size);
		
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
				
			}
			
			if (sigma_node != 0)
				local_S *= 0.5*t_step;

			jj=kk=0;
			// if (nid==0)
				// std::cout << std::endl << local_E << std::endl << std::endl;
			/* Handling pec edges */
			auto local_E_nondirichlet = local_E;
			for (auto j = global_i.begin(); j != global_i.end(); j++)
			{	
				if (classify_edges[*j] == 4)
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
				// else if (classify_edges[*j] == 2)
				// {
					
				// }
				jj++;
			}
			
			// std::cout << std::endl << local_Id << std::endl << std::endl;
			
			// auto local_H = (local_E+local_S).inverse(); 
			Eigen::MatrixXd local_H = (local_E+local_S).llt().solve(local_Id); //faster than inverse
			Eigen::MatrixXd local_E_inv = (local_E_nondirichlet).llt().solve(local_Id);
			auto local_P = local_H*(local_E - local_S);
	
			// std::cout << std::endl << local_P << std::endl << std::endl;
			// std::cout << std::endl << local_H << std::endl << std::endl;
			
			uint32_t offset;
			bool is_frac=false;
			if (dual_is_fractured[nid]>0)
			{
				if (dual_is_fractured[nid] >= Q_size+1)
					std::cout << nid << "---" << dual_is_fractured[nid] << "---" << global_i.size() << std::endl;
				is_frac=true;
				offset = dual_is_fractured[nid]-1;
			}
			// else
				// std::cout << std::endl << local_P << std::endl << std::endl;
			jj=kk=0;
			
			for (auto j = global_i.begin(); j != global_i.end(); j++)
			{
				switch (classify_edges[*j])
				{
					case 1 :
					{
						for (auto k = global_i.begin(); k != global_i.end(); k++)
						{
							if (classify_edges[*k] !=4 && local_H(jj,kk)!=0)
								H_trip.push_back(double_triplet(h_index[*j],*k,local_H(jj,kk)));
							kk++;
						}
						
						M_trip.push_back(double_triplet(*j,h_index[*j],double(1)));
						break;
					}
					case 2 :
					{
						P_p[p_index[*j]]=local_P(jj,jj);
						for (auto k = global_i.begin(); k != global_i.end(); k++)
						{
							/*if (local_P(jj,kk)!=0)
							{
								if (classify_edges[*k]==1)
								{
									P_trip.push_back(double_triplet(p_index[*j],h_index[*k],local_P(jj,kk)));
									std::cout << "I didn't think this would ever happen... Go figure!" << std::endl;
								}
								else if (classify_edges[*k]==2)
									P_trip.push_back(double_triplet(p_index[*j],H_size+p_index[*k],0.5*local_P(jj,kk)));
								else if (classify_edges[*k]==3)
								{
									P_trip.push_back(double_triplet(p_index[*j],H_size+P_size+offset+kk,local_P(jj,kk)));
									std::cout << "I didn't think this would ever happen... Go figure!" << std::endl;
								}
							}*/

							if (local_H(jj,kk)!=0)
								Mp_trip.push_back(double_triplet(p_index[*j],*k,local_H(jj,kk)));
							
							kk++;
						}
						
						M_trip.push_back(double_triplet(*j,H_size+p_index[*j],double(1)));
						break;
					}
					case 3 :
					{
						M_trip.push_back(double_triplet(*j,H_size+P_size+offset+jj,double(1)));
						
						break;
					}
					case 4 :
					{
						M_trip.push_back(double_triplet(*j,H_size+P_size+Q_size+boundary_index[*j],double(1)));
						break;
					}
					default :
					{
						std::cout << "Edge domain class = " << classify_edges[*j] << ". This is unexpected!" << std::endl;
						break;
					}
				}
				// }
				
				if (is_frac)
				{
					if (classify_edges[*j] != 4)
					{
						for (kk = 0; kk < local_E_size; kk++)
						{
							if (classify_edges[global_i[kk]] != 4)
							{
								if (local_P(jj,kk) != 0)
									Q_trip.push_back(double_triplet(offset+jj,offset+kk,local_P(jj,kk)));
							
								if (local_H(jj,kk) != 0)
								{
									Mq_trip.push_back(double_triplet(offset+jj,global_i[kk],local_H(jj,kk)));
									if (offset+jj > Q_size-1 || global_i[kk] >= edges_size())
										std::cout << Q_size << " " << global_i.size() << " " << edges_size() << " " << global_i[kk] << std::endl;
								}
							}
						}						
					}
				}
				
				for (uint32_t ll=0; ll < local_E_size; ll++)
					if (local_E_inv(jj,ll) != 0)
						Einv_trip.push_back(double_triplet(global_i[jj],global_i[ll],local_E_inv(jj,ll)));
				
				jj++;
				kk=0;
			}
		}
		
		
		// t_material.toc();
		// std::cout << "done - " << t_material << " seconds" << std::endl;
		
		// t_material.tic();
		U_frac_size = H_size + P_size + Q_size + B_size;
		F_frac_size = N_size + R_size + S_size;
		
		this->P_p=std::move(P_p);
		
		Eigen::SparseMatrix<double> H(H_size,edges_size()), N(N_size,surfaces_size()), Einv(edges_size(),edges_size());
		Eigen::SparseMatrix<double> M(edges_size(),U_frac_size), Mq(Q_size,edges_size()), Mp(P_size,edges_size());
		Eigen::SparseMatrix<double> T(surfaces_size(),F_frac_size),Tr(R_size,surfaces_size()),Ts(S_size,surfaces_size());
		Eigen::SparseMatrix<double> /*P(P_size,U_frac_size),*/Q(Q_size,Q_size);
		Eigen::SparseMatrix<double> R(R_size,F_frac_size),S(S_size,S_size);
			
		add_to_sparse ass;
		overwrite_to_sparse oss;
		
		T.setFromTriplets(T_trip.begin(),T_trip.end(), oss);
		// std::cout << "1" << std::endl;
		M.setFromTriplets(M_trip.begin(),M_trip.end(), oss);
		// std::cout << "2" << std::endl;
		N.setFromTriplets(N_trip.begin(),N_trip.end(), ass);
		// std::cout << "3" << std::endl;
		H.setFromTriplets(H_trip.begin(),H_trip.end(), ass);
		// std::cout << "4" << std::endl;
		Mp.setFromTriplets(Mp_trip.begin(),Mp_trip.end(), ass);
		// std::cout << "5" << std::endl;
		Mq.setFromTriplets(Mq_trip.begin(),Mq_trip.end(), ass);
		// std::cout << "6" << std::endl;
		N.setFromTriplets(N_trip.begin(),N_trip.end(), ass);
		// std::cout << "7" << std::endl;
		P.setFromTriplets(P_trip.begin(),P_trip.end(), ass);
		// std::cout << "8" << std::endl;
		Q.setFromTriplets(Q_trip.begin(),Q_trip.end(), ass);
		// std::cout << "9" << std::endl;
		R.setFromTriplets(R_trip.begin(),R_trip.end(), ass);
		// std::cout << "10" << std::endl;
		S.setFromTriplets(S_trip.begin(),S_trip.end(), ass);
		// std::cout << "11" << std::endl;
		Tr.setFromTriplets(Tr_trip.begin(),Tr_trip.end(), ass);
		// std::cout << "12" << std::endl;
		Ts.setFromTriplets(Ts_trip.begin(),Ts_trip.end(), ass);
		// std::cout << "13" << std::endl;
		Einv.setFromTriplets(Einv_trip.begin(),Einv_trip.end(), ass);
		// std::cout << "14" << std::endl;
		// Eigen::MatrixXd dMat = Eigen::MatrixXd(H);
		// std::ofstream debug_file("new_massmatrix.dat");
		// debug_file << dMat << std::endl;
		// debug_file.close();
		
		this->H=std::move(H); this->Einv=std::move(Einv);
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
		return (pts[std::get<0>(volumes[v])]+pts[std::get<1>(volumes[v])]+pts[std::get<2>(volumes[v])]+pts[std::get<3>(volumes[v])])/4;
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
	uint32_t									input_line, H_size, Q_size, P_size, B_size, N_size, R_size, S_size, U_frac_size, F_frac_size;
	Eigen::SparseMatrix<double> 				C,H,M,Mq,Mp,N,P,Q,R,S,T,Tr,Ts,Einv;
	Eigen::VectorXd								U,F,I,P_p;
	std::mutex									meshlock;
	std::map<uint32_t,Output>					Outputs;
	std::map<uint32_t,Simulation>				Simulations;
	std::map<uint32_t,Source>					Sources;								/* a std::map works because every time I use the [] operator on an undefined material */
	std::map<uint32_t,Material>					Materials; 								/* (or source), the default constructor makes it empty space (or null source) */
	std::map<uint32_t,BoundaryCondition>		BCs;
	std::map<uint32_t,Mesh>						Meshes;
	std::vector<double>                         CellVolumes, probe_numeric_times;
	std::vector<std::vector<double>>			probe_numeric_xvalues, probe_numeric_yvalues, probe_numeric_zvalues, analytic_value_vector;
	// std::vector<bool>							is_dirichlet;
	std::vector<uint32_t> 						vol_material, boundary_index, h_index, p_index, n_index, r_index, edge_bids;
	std::vector<uint32_t>				        edge_bcs, face_bcs, probe_elem;
	std::vector<uint8_t>						classify_edges, classify_surfaces;
	std::vector<std::vector<uint32_t>>          associated_volumes, dual_star_offsets;
	std::vector<std::vector<uint32_t>>			associated_frac_edges, associated_p_edges, associated_h_edges, associated_bnd_edges; 	
	std::vector<std::vector<uint32_t>>			edge_src, face_src;
	std::vector<bool> 							bnd_nodes;
	std::vector<volume_type> 					volumes;
	std::vector<surface_type> 					surfaces;
	std::vector<uint32_t>                       dual_is_fractured, primal_is_fractured;
	std::vector<edge_type> 						edges;
	std::vector<Eigen::Vector3d>	 			pts, edge_bars, face_bars, probe_points;
	std::vector<cluster_list>    				nte_list, etn_list, etf_list, fte_list, ftv_list, vtf_list;								
	double                                      t_step, min_h, average_diameter;
	uint32_t									loaded_mesh_label, current_simulation;
	DBfile 										*_siloDb=NULL;
	// std::array<std::vector<uint32_t>,20>		sources_by_label; 						/* Each label has a vector containing all the sources active on that label */
};
