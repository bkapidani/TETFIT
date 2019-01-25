// file Discretization.hpp
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
#include "Utilities.hpp" //contains also all includes to c++ std libraries, BLAZE and SILOs
#include "Source.hpp"
#include "BoundaryCondition.hpp"
#include "Material.hpp"
#include "Mesh.hpp"
#include "Solid.hpp"
#include "Simulation.hpp"
#include "Output.hpp"
//~ #include "ConjugateGradientSolver.hpp"

class Discretization
{   
   public:
   Discretization(std::string inputfile)
   {
      std::string line;
      input_line = 1;
      root = 0;
      std::ifstream ReadFile;//(inputfile.c_str());
      ReadFile.open(inputfile.c_str());
      bool in_definition = false;
      // bool freq_not_set = true;
      // bool mesh_added = false;
      
      
      char action[10], token[20], value[64];
      std::string thing_being_defined;
      uint32_t definition_label,input_line;

      std::vector<uint32_t> BndToBC(100,0);   
      // std::cout << "ciao!" << std::endl;
      input_line = 1;
      while(getline(ReadFile,line))
      {
         // std::cout << "ciao!" << std::endl;
         ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;
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
                        definition_label = std::stoi(val);
                        in_definition = true;
                        
                        //~ std::cout << definition_label << std::endl;
                        
                        if (thing_being_defined == "material")
                        {
                           Material matt;
                           Materials.insert(std::pair<uint32_t,Material>(definition_label,matt));
                        }
                        else if (thing_being_defined == "bc")
                        {
                           BoundaryCondition bici;
                           BCs.insert(std::pair<uint32_t,BoundaryCondition>(definition_label,bici));
                        }
                        else if (thing_being_defined == "mesh")
                        {
                           Mesh mascia;
                           Meshes.insert(std::pair<uint32_t,Mesh>(definition_label,mascia));
                        }
                        else if (thing_being_defined == "refinement")
                        {
                           Refinement raffa;
                           Refinements.insert(std::pair<uint32_t,Refinement>(definition_label,raffa));
                        }
                        else if (thing_being_defined == "source")
                        {
                           Source saorsie;
                           Sources.insert(std::pair<uint32_t,Source>(definition_label,saorsie));
                        }
                        else if (thing_being_defined == "simulation")
                        {
                           Simulation simo;
                           Simulations.insert(std::pair<uint32_t,Simulation>(definition_label,simo));
                        }
                        else if (thing_being_defined == "output")
                        {
                           Output pat;
                           Outputs.insert(std::pair<uint32_t,Output>(definition_label,pat));
                        }
                     }
                           ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

                  }
                        ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

               }
               else if (instr == "END")
               {
                  if (tok != thing_being_defined || (static_cast<uint32_t>(std::atoi(value)) != definition_label))
                     MyThrow(input_line,end_wo_define);
                  else 
                     in_definition = false;
                  
                  //~ std::cout << thing_being_defined << std::endl;
               // std::cout << "did not crash at script line " << input_line << std::endl;
                  
                  if (thing_being_defined == "bc")
                     for (auto surfz : BCs[definition_label].Surfaces())
                        BndToBC[uint32_t(surfz)] = definition_label;
                  
               
               // std::cout << "again did not crash at script line " << input_line << std::endl;
               }
               else if (instr != "SET")
                  MyThrow(input_line,unknown_instruction);
               else
               {
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
                  else if (thing_being_defined == "refinement")
                     Refinements[definition_label].SetParam(input_line,tok,val);
                  else if (thing_being_defined == "source")
                  {
                     ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

                     Sources[definition_label].SetParam(input_line,tok,val);
                  }
                  else if (thing_being_defined == "simulation")
                  {
                     Simulations[definition_label].SetParam(input_line,tok,val);
                     if (tok == "method")
                        method_line = input_line;
                  }
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
      
      // for (auto bbc : BCs)
         // std::cout << bbc.first << " " << bbc.second.Type() << std::endl;
      // std::cout << std::endl;
      
      this->BndToBC = std::move(BndToBC);
      ReadFile.close();
      
      ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

      // Run();
   }
   
   Discretization(const Discretization& disc)
   {
      root = 0;
      input_line = 1;
      this->method_line = disc.method_line;
      
      ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

      this->Materials    = disc.Materials;
      this->BCs          = disc.BCs;
      this->Meshes       = disc.Meshes;
      
      ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

      this->Sources       = disc.Sources;
      
      ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

      this->Simulations    = disc.Simulations;
      this->Outputs       = disc.Outputs;
      this->Solids       = disc.Solids;
      this->Refinements   = disc.Refinements;
      this->BndToBC      = disc.BndToBC;
      
      ////std::cout << __FILE__ << ":" << __LINE__ << std::endl;

      // for (auto bbc : BCs)
         // std::cout << bbc.first << " " << bbc.second.Type() << std::endl;
      // std::cout << std::endl;
   }

   //~ void Run(void)
   //~ {
      //~ timecounter t_sim;
      //~ t_sim.tic();
      
      //~ for (auto sims : Simulations)
         //~ RunSimulation(sims.second, sims.first);
      //~ t_sim.toc();
      
      //~ std::cout << std::endl;
      //~ std::cout << "----------- All simulations ran. Elapsed time: " << t_sim << " seconds ----------" << std::endl;
   //~ }
   
   void RunSimulation(const Simulation& s, const uint32_t& sim_label)
   {
      std::cout  << std::endl << std::endl;
      
      DateAndTime();
      
      std::cout <<"------------------------ Running Time Domain Simulation ------------------------" << std::endl << std::endl;
      
      bool probes_out_of_mesh = false;
      //~ bool dipoles_out_of_mesh = false;
      //~ double current_time;
      this->excitation_freq = Sources[*(Simulations[current_simulation].Src().begin())].GetFreq();
      auto m = Meshes[s.MeshLabel()];
      auto meth = s.Method();
      auto o = &Outputs[s.Output()];
      auto mod_out = (*o).Mode();
      timecounter t_preproc, step_cost;
      std::vector<double> Losses;


      std::cout << "Preprocessing... ";
      std::cout.flush();
      
      t_preproc.tic();

      possibly_unstable = false;
      current_simulation = sim_label;
      (*o).Initialize();
      simulation_time = s.Time();
      
      /***********************************MESH PARSING************************************/
      
      ReadMesh(m);
      
      /***********************************************************************************/
            
      double c_0 = 1; //vacuum speed of light in natural units
      t_step = 0.5*min_height/c_0;
      
      if (s.Courant()>1)
         possibly_unstable = true;
      else
         t_step *= s.Courant();
      
      // time counters
      double export_time_average,step_time_average,bcs_time_average,mag_time_average,ele_time_average;
      export_time_average=bcs_time_average=mag_time_average=ele_time_average=step_time_average=0;
      
      t_preproc.toc();
      std::cout << " done (" << t_preproc << " seconds)" << std::endl;
      
      std::cout << "Preparing output... ";
      std::cout.flush();
      t_preproc.tic();
      
      uint32_t NP = 0;
      std::vector<std::vector<std::array<double,3UL>>> probe_numeric_E,probe_numeric_H;
      std::vector<std::vector<std::array<double,3UL>>> probe_analytic_E,probe_analytic_H;
      
      if ( mod_out == "probefield" )
      {
         //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         
         NP=(*o).Nprobes();
         probe_elem.reserve(NP);
         
         for (uint32_t p=0; p<NP; ++p)
         {
            auto pv = (*o).Probe(p);
            auto elem_index = FindProbe(pv);
            
            if (elem_index.first < elements.size())
            {
               probe_elem.push_back(elem_index);
               probe_index.push_back(p);
            }else 
               probes_out_of_mesh = true;
         }
         std::ofstream op((*o).Name().c_str());
         op.close();
      }
      
      t_preproc.toc();
      std::cout << " done (" << t_preproc << " seconds)" << std::endl;

      if (probes_out_of_mesh)
         std::cout << "WARNING: one or more field probes were set out of the mesh! (And will be therefore ignored)" << std::endl;
      if (possibly_unstable)
         std::cout << "WARNING: Courant factor > 1 leads to instability! (automatically set to 1)" << std::endl;
      
      std::cout << std::endl << "Simulation parameters:"       << std::endl;
      std::cout << std::setw(20) << "Mesh: "                   << std::setw(20) << m.FileName() << std::endl;
      std::cout << std::setw(20) << "Domain bounding box: "    << std::setw(10) << "{" << xmin << "," << xmax << "} X "
                                                               << "{" << ymin << "," << ymax << "} X "
                                                               << "{" << zmin << "," << zmax << "}" << std::endl;
      std::cout << std::setw(20) << "Simulation time: "
                << std::setw(20) << simulation_time << " sec" << std::endl;
      std::cout << std::setw(20) << "Time step: "
                << std::setw(20) << t_step << " sec" << std::endl;
      std::cout << std::setw(20) << "Elements: "
                << std::setw(20) << elements.size() << std::endl;
      std::cout << std::setw(20) << "Points:   "
                << std::setw(20) << pts.size() << std::endl;
      std::cout << std::setw(20) << "Unknowns: "
                << std::setw(20) << facets.size()+edges.size() << std::endl;

      uint32_t i=0;
      blaze::StaticMatrix<double,4UL,6UL> C;
      blaze::StaticMatrix<double,6UL,4UL> CT;
      
      C  = {{ 0, 0, 1, 0,-1, 1},
            { 0, 1, 0,-1, 0, 1},
            { 1, 0, 0,-1, 1, 0},
            { 1,-1, 1, 0, 0, 0}};
            
      CT = {{ 0, 0, 1, 1},
            { 0, 1, 0,-1},
            { 1, 0, 0, 1},
            { 0,-1,-1, 0},
            {-1, 0, 1, 0},
            { 1, 1, 0, 0}};
      
      // 2 forms (flux global quantities) for B and D fields
      BlazeDVec B2_form(edges.size()),H1_form(edges.size());
      B2_form = H1_form = 0;
      
      timecounter step_counter, tdbg;
      double current_time=0;
      double peak=0.5;
      uint32_t N_of_steps = uint32_t(simulation_time/t_step);
      
      i=0;
      bool init_conds=true;
         // Ampere-Maxwell:
      if (init_conds)
      {
         BlazeDVec new_B2form(edges.size());
         new_B2form = 0;

         for (uint32_t vol=0;vol<elements.size();++vol)
         {
            for (uint32_t ee=0; ee<6; ++ee)
            {
               auto edgar = elem_edges[vol][ee];
               auto n0 = std::get<0>(edges[edgar]);
               auto n1 = std::get<1>(edges[edgar]);
               auto p0 = pts[n0];
               auto p1 = pts[n1];
               
               std::array<double,3UL> mp0({(0.75*p0[0]+0.25*p1[0]),
                                           (0.75*p0[1]+0.25*p1[1]),
                                           (0.75*p0[2]+0.25*p1[2])});
                                           
               std::array<double,3UL> mp1({(0.25*p0[0]+0.75*p1[0]),
                                           (0.25*p0[1]+0.75*p1[1]),
                                           (0.25*p0[2]+0.75*p1[2])});
               
               Blaze3Vec v0({0.5*(p1[0]-p0[0]),0.5*(p1[1]-p0[1]),0.5*(p1[2]-p0[2])});
               
               double sv0 = std::exp(-100*(std::pow(mp0[1]-peak,2)+std::pow(mp0[2]-peak,2)));
               double sv1 = std::exp(-100*(std::pow(mp1[1]-peak,2)+std::pow(mp1[2]-peak,2)));
               
               
               Blaze3Vec fv0({sv0,0,0});
               Blaze3Vec fv1({sv1,0,0});
               
               
               uint32_t i0 = dual_index[n0]+global_node_offset_pairs[edgar][0];
               uint32_t i1 = dual_index[n1]+global_node_offset_pairs[edgar][1];
               local_F[i0] = (fv0,v0);
               local_F[i1] = (fv1,v0);
               H1_form[edgar] = local_F[i0]+local_F[i1];
            }
            
            auto LU = blaze::subvector(local_U,4*vol,4);
            LU = 0;
                        
            auto p0 = std::get<0>(elements[vol]);
            auto p1 = std::get<1>(elements[vol]);
            auto p2 = std::get<2>(elements[vol]);
            auto p3 = std::get<3>(elements[vol]);
            
            Blaze3Mat primmat;
            primmat = { {pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1], pts[p1][2]-pts[p0][2]},
                        {pts[p2][0]-pts[p0][0], pts[p2][1]-pts[p0][1], pts[p2][2]-pts[p0][2]},
                        {pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1], pts[p3][2]-pts[p0][2]} };
                        
            auto volume = det(primmat);
            std::array<double,4> coeffs({1,-1,1,-1});
            if (volume<0)
            {
               coeffs[0]=-coeffs[0];
               coeffs[1]=-coeffs[1];
               coeffs[2]=-coeffs[2];
               coeffs[3]=-coeffs[3];
            }
            
            std::array<double,3UL> vb0({0.25*(pts[p0][0]+pts[p1][0]+pts[p2][0]+pts[p3][0]),
                                        0.25*(pts[p0][1]+pts[p1][1]+pts[p2][1]+pts[p3][1]),
                                        0.25*(pts[p0][2]+pts[p1][2]+pts[p2][2]+pts[p3][2])});
                  
            std::array<double,3UL> fb0({(pts[p1][0]+pts[p2][0]+pts[p3][0])/double(3),
                                        (pts[p1][1]+pts[p2][1]+pts[p3][1])/double(3),
                                        (pts[p1][2]+pts[p2][2]+pts[p3][2])/double(3)});
            std::array<double,3UL> fb1({(pts[p0][0]+pts[p2][0]+pts[p3][0])/double(3),
                                        (pts[p0][1]+pts[p2][1]+pts[p3][1])/double(3),
                                        (pts[p0][2]+pts[p2][2]+pts[p3][2])/double(3)});
            std::array<double,3UL> fb2({(pts[p0][0]+pts[p1][0]+pts[p3][0])/double(3),
                                        (pts[p0][1]+pts[p1][1]+pts[p3][1])/double(3),
                                        (pts[p0][2]+pts[p1][2]+pts[p3][2])/double(3)});
            std::array<double,3UL> fb3({(pts[p0][0]+pts[p1][0]+pts[p2][0])/double(3),
                                        (pts[p0][1]+pts[p1][1]+pts[p2][1])/double(3),
                                        (pts[p0][2]+pts[p1][2]+pts[p2][2])/double(3)});
            
            std::array<double,3UL> mp0({(vb0[0]+fb0[0])*0.5,(vb0[1]+fb0[1])*0.5,(vb0[2]+fb0[2])*0.5});
            std::array<double,3UL> mp1({(vb0[0]+fb1[0])*0.5,(vb0[1]+fb1[1])*0.5,(vb0[2]+fb1[2])*0.5});
            std::array<double,3UL> mp2({(vb0[0]+fb2[0])*0.5,(vb0[1]+fb2[1])*0.5,(vb0[2]+fb2[2])*0.5});
            std::array<double,3UL> mp3({(vb0[0]+fb3[0])*0.5,(vb0[1]+fb3[1])*0.5,(vb0[2]+fb3[2])*0.5});
            
            Blaze3Vec dv0({-(vb0[0]-fb0[0]),-(vb0[1]-fb0[1]),-(vb0[2]-fb0[2])});
            Blaze3Vec dv1({-(vb0[0]-fb1[0]),-(vb0[1]-fb1[1]),-(vb0[2]-fb1[2])});
            Blaze3Vec dv2({-(vb0[0]-fb2[0]),-(vb0[1]-fb2[1]),-(vb0[2]-fb2[2])});
            Blaze3Vec dv3({-(vb0[0]-fb3[0]),-(vb0[1]-fb3[1]),-(vb0[2]-fb3[2])});
                        
            double sv0 = std::exp(-100*(/*std::pow(mp0[0]-peak,2)+*/std::pow(mp0[1]-peak,2)+std::pow(mp0[2]-peak,2)));
            double sv1 = std::exp(-100*(/*std::pow(mp1[0]-peak,2)+*/std::pow(mp1[1]-peak,2)+std::pow(mp1[2]-peak,2)));
            double sv2 = std::exp(-100*(/*std::pow(mp2[0]-peak,2)+*/std::pow(mp2[1]-peak,2)+std::pow(mp2[2]-peak,2)));
            double sv3 = std::exp(-100*(/*std::pow(mp3[0]-peak,2)+*/std::pow(mp3[1]-peak,2)+std::pow(mp3[2]-peak,2)));
            
            Blaze3Vec fv0({sv0,0,0});
            Blaze3Vec fv1({sv1,0,0});
            Blaze3Vec fv2({sv2,0,0});
            Blaze3Vec fv3({sv3,0,0});
            
            LU[0] = coeffs[0]*(fv0,dv0);
            LU[1] = coeffs[1]*(fv1,dv1);
            LU[2] = coeffs[2]*(fv2,dv2);
            LU[3] = coeffs[3]*(fv3,dv3);
         
            blaze::StaticVector<double,6UL> local_Phi;
            local_Phi = CT * LU;
            for (uint8_t k=0;k<6;++k)
               new_B2form[elem_edges[vol][k]] = new_B2form[elem_edges[vol][k]] + 
                                                local_Phi[k];
         }
         
         B2_form = std::move(new_B2form);
      }
      
      while (current_time < simulation_time)
      {
         step_cost.tic();
         tdbg.tic();
         
         if ((*o).AllowPrint(current_time/*-t_step*/))
            PlotFields(mod_out, current_time/*-t_step*/,uint32_t(i));

         tdbg.toc();
         export_time_average += (duration_cast<duration<double>>(tdbg.elapsed())).count();
         //~ tdbg.tic();
         
         //sources
         //~ for (auto be : bnd_primal)
            //~ B2_form[be] += ComputeEfieldSource(be,current_time);
         //~ tdbg.toc();
         //~ bcs_time_average += (duration_cast<duration<double>>(tdbg.elapsed())).count();
         tdbg.tic();
         
         // Faraday:
         BlazeDVec new_H1form(edges.size());
         new_H1form = 0;
         for (uint32_t p=0;p<pts.size();++p)
         {
            BlazeDVec local_Phi(n_e[p].size());
            auto LF = blaze::subvector(local_F,dual_index[p],n_e[p].size());
            
            for (uint8_t k=0;k<n_e[p].size();++k)
            {
               //~ auto p0 = std::get<0>(edges[n_e[p][k]]);
               local_Phi[k] = B2_form[n_e[p][k]];
               //~ /*local_Phi[k] = (p0 == p) ?  B2_form[n_e[p][k]] :
                                          //~ -B2_form[n_e[p][k]] ;*/
            }
            
            LF -= t_step*(NU[p]*local_Phi);
            
            for (uint8_t k=0;k<n_e[p].size();++k)
            {
               //~ auto p0 = std::get<0>(edges[n_e[p][k]]);
               new_H1form[n_e[p][k]] += LF[k];
               //~ new_H1form[n_e[p][k]] += (p0 == p) ? LF[k] : -LF[k];
            }
         }
         
         H1_form = std::move(new_H1form);
         
         tdbg.toc();
         mag_time_average += (duration_cast<duration<double>>(tdbg.elapsed())).count();
         tdbg.tic();
         
         // Ampere-Maxwell:
         BlazeDVec new_B2form(edges.size());
         new_B2form = 0;
         for (uint32_t v=0; v<elements.size(); ++v)
         {
            blaze::StaticVector<double,6UL> local_F0;
            
            for (uint8_t k=0;k<6;++k)
            {
               local_F0[k] = H1_form[elem_edges[v+0][k]];
            }
            
            auto local_Eta0 = blaze::submatrix(ETA,4*(v+0),0UL,4UL,4UL);
            auto LU0 = blaze::subvector(local_U,4*(v+0),4);
            
            LU0 += blaze::eval(t_step*(local_Eta0*(C * local_F0)));
            blaze::StaticVector<double,6UL> local_Phi0 = blaze::eval(CT * LU0);

            for (uint8_t k=0;k<6;++k)
            {
               new_B2form[elem_edges[v+0][k]] += local_Phi0[k];
            }
         }
         
         B2_form = std::move(new_B2form);
         
         tdbg.toc();
         //~ std::cout << "profiling avx: " << tdbg << std::endl;
         ele_time_average += (duration_cast<duration<double>>(tdbg.elapsed())).count();
         
         step_cost.toc();
         step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
         
         if (i !=0 && i % 137 == 0) //arbitrary, just a nod at the fine structure constant
            std::cout << "-----------" << "Progress: " << std::setw(2) << 100*i/N_of_steps 
                      << "% done in " << std::setw(9) << step_time_average << "s, " 
                      << std::setw(8) << step_time_average/double(i) << std::setw(7) << " s/step" 
                      << "-----------" << std::endl;
         ++i;
         current_time = double(i)*t_step;
      }
      
      //~ std::stringstream fname;
      //~ fname << (*o).Name() << ".png";
      //~ plt::save(fname.str());
      
      std::cout << std::endl << "Simulation statistics:"       << std::endl;
      std::cout << std::setw(20) << "Average step cost in seconds:           "
                << std::setw(20) << step_time_average/double(i) << std::endl;
      std::cout << std::setw(20) << "Average plot cost in seconds:           "
                << std::setw(20) << export_time_average/double(i) << std::endl;
      std::cout << std::setw(20) << "Average Faraday law cost in seconds:    "
                << std::setw(20) << mag_time_average/double(i) << std::endl;
      std::cout << std::setw(20) << "Average Ampere-Maxwell cost in seconds: "
                << std::setw(20) << ele_time_average/double(i) << std::endl;
      return;
   }
   
   void PlotFields(const std::string s, double t, uint32_t i)
   {
      // Eigen::Vector3d num_val;
      if (s == "probefield")
      {
         auto op = &Outputs[Simulations[current_simulation].Output()];
         std::array<std::vector<double>,6> ordinata;
         if (t<=0)
         {
            if ((*op).XG())
            {
               for (uint32_t p=0; p<probe_index.size(); ++p)
               {
                  auto ppt = (*op).Probe(probe_index[p]);
                  abscissa.push_back(ppt[0]);
               }
            }
            else if ((*op).YG())
            {
               for (uint32_t p=0; p<probe_index.size(); ++p)
               {
                  auto ppt = (*op).Probe(probe_index[p]);
                  abscissa.push_back(ppt[1]);
               }
            }
            else if ((*op).ZG())
            {
               for (uint32_t p=0; p<probe_index.size(); ++p)
               {
                  auto ppt = (*op).Probe(probe_index[p]);
                  abscissa.push_back(ppt[2]);
               }
            }
            else
            {
               uint32_t p = probe_index.size();
               
               std::vector<std::vector<double>> exs(p);
               std::vector<std::vector<double>> eys(p);
               std::vector<std::vector<double>> ezs(p);
               std::vector<std::vector<double>> hxs(p);
               std::vector<std::vector<double>> hys(p);
               std::vector<std::vector<double>> hzs(p);
               
               this->exs = std::move(exs);
               this->eys = std::move(eys);
               this->ezs = std::move(ezs);
               this->hxs = std::move(hxs);
               this->hys = std::move(hys);
               this->hzs = std::move(hzs);
            }
         }
         
         
         //~ std::ofstream probe_file((*op).Name().c_str(),std::ofstream::out | std::ofstream::app);
         for (uint32_t p=0; p<probe_index.size(); ++p)
         {
            auto vol = probe_elem[p].first;
            Blaze3Vec num_ele(0), num_mag(0);
            
            auto p0 = std::get<0>(elements[vol]);
            auto p1 = std::get<1>(elements[vol]);
            auto p2 = std::get<2>(elements[vol]);
            auto p3 = std::get<3>(elements[vol]);
            
            std::array<double,3UL> vb0({0.25*(pts[p0][0]+pts[p1][0]+pts[p2][0]+pts[p3][0]),
                                        0.25*(pts[p0][1]+pts[p1][1]+pts[p2][1]+pts[p3][1]),
                                        0.25*(pts[p0][2]+pts[p1][2]+pts[p2][2]+pts[p3][2])});
                  
            std::array<double,3UL> fb0({(pts[p1][0]+pts[p2][0]+pts[p3][0])/double(3),
                                        (pts[p1][1]+pts[p2][1]+pts[p3][1])/double(3),
                                        (pts[p1][2]+pts[p2][2]+pts[p3][2])/double(3)});
            std::array<double,3UL> fb1({(pts[p0][0]+pts[p2][0]+pts[p3][0])/double(3),
                                        (pts[p0][1]+pts[p2][1]+pts[p3][1])/double(3),
                                        (pts[p0][2]+pts[p2][2]+pts[p3][2])/double(3)});
            std::array<double,3UL> fb2({(pts[p0][0]+pts[p1][0]+pts[p3][0])/double(3),
                                        (pts[p0][1]+pts[p1][1]+pts[p3][1])/double(3),
                                        (pts[p0][2]+pts[p1][2]+pts[p3][2])/double(3)});
            std::array<double,3UL> fb3({(pts[p0][0]+pts[p1][0]+pts[p2][0])/double(3),
                                        (pts[p0][1]+pts[p1][1]+pts[p2][1])/double(3),
                                        (pts[p0][2]+pts[p1][2]+pts[p2][2])/double(3)});
            
            //~ std::array<double,3UL> eb0({0.5*(pts[p0][0]+pts[p1][0]),
                                        //~ 0.5*(pts[p0][1]+pts[p1][1]),
                                        //~ 0.5*(pts[p0][2]+pts[p1][2])});
            //~ std::array<double,3UL> eb1({0.5*(pts[p0][0]+pts[p2][0]),
                                        //~ 0.5*(pts[p0][1]+pts[p2][1]),
                                        //~ 0.5*(pts[p0][2]+pts[p2][2])});
            //~ std::array<double,3UL> eb2({0.5*(pts[p1][0]+pts[p2][0]),
                                        //~ 0.5*(pts[p1][1]+pts[p2][1]),
                                        //~ 0.5*(pts[p1][2]+pts[p2][2])});
            //~ std::array<double,3UL> eb3({0.5*(pts[p0][0]+pts[p3][0]),
                                        //~ 0.5*(pts[p0][1]+pts[p3][1]),
                                        //~ 0.5*(pts[p0][2]+pts[p3][2])});
            //~ std::array<double,3UL> eb4({0.5*(pts[p1][0]+pts[p3][0]),
                                        //~ 0.5*(pts[p1][1]+pts[p3][1]),
                                        //~ 0.5*(pts[p1][2]+pts[p3][2])});
            //~ std::array<double,3UL> eb5({0.5*(pts[p2][0]+pts[p3][0]),
                                        //~ 0.5*(pts[p2][1]+pts[p3][1]),
                                        //~ 0.5*(pts[p2][2]+pts[p3][2])});
            
            //~ auto ae0 = 0.5*((vb0-eb0)%(fb2-fb3));
            //~ auto ae1 = 0.5*((vb0-eb1)%(fb3-fb1));
            //~ auto ae2 = 0.5*((vb0-eb2)%(fb0-fb3));
            //~ auto ae3 = 0.5*((vb0-eb3)%(fb1-fb2));
            //~ auto ae4 = 0.5*((vb0-eb4)%(fb2-fb0));
            //~ auto ae5 = 0.5*((vb0-eb5)%(fb0-fb1));
            
            
            
            Blaze3Mat primmat,dualmat;
            primmat = { {pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1], pts[p1][2]-pts[p0][2]},
                        {pts[p2][0]-pts[p0][0], pts[p2][1]-pts[p0][1], pts[p2][2]-pts[p0][2]},
                        {pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1], pts[p3][2]-pts[p0][2]} };
                        
            auto volume = det(primmat);
            std::array<double,4> coeffs({1,-1,1,-1});
            if (volume<0)
            {
               coeffs[0]=-coeffs[0];
               coeffs[1]=-coeffs[1];
               coeffs[2]=-coeffs[2];
               coeffs[3]=-coeffs[3];
            }
            std::array<double,3> dual_inds,prim_inds;
            
            switch(probe_elem[p].second)
            {
               case 0 :
               {
                  prim_inds[0] = 4*vol+1; prim_inds[1] = 4*vol+2; prim_inds[2] = 4*vol+3;
                  dual_inds[0] = dual_index[p0]+global_node_offset_pairs[elem_edges[vol][0]][0];
                  dual_inds[1] = dual_index[p0]+global_node_offset_pairs[elem_edges[vol][1]][0];
                  dual_inds[2] = dual_index[p0]+global_node_offset_pairs[elem_edges[vol][3]][0];
                  
                  //~ primmat = { {pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1], pts[p1][2]-pts[p0][2]},
                              //~ {pts[p2][0]-pts[p0][0], pts[p2][1]-pts[p0][1], pts[p2][2]-pts[p0][2]},
                              //~ {pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1], pts[p3][2]-pts[p0][2]} };
                  
                  dualmat = { {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                              {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])},
                              {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                  break;
               }
               case 1 :
               {
                  prim_inds[0] = 4*vol+0; prim_inds[1] = 4*vol+2; prim_inds[2] = 4*vol+3;
                  dual_inds[0] = dual_index[p1]+global_node_offset_pairs[elem_edges[vol][0]][1];
                  dual_inds[1] = dual_index[p1]+global_node_offset_pairs[elem_edges[vol][2]][0];
                  dual_inds[2] = dual_index[p1]+global_node_offset_pairs[elem_edges[vol][4]][0];
                  
                  primmat = { {pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1], pts[p1][2]-pts[p0][2]},
                              {pts[p2][0]-pts[p1][0], pts[p2][1]-pts[p1][1], pts[p2][2]-pts[p1][2]},
                              {pts[p3][0]-pts[p1][0], pts[p3][1]-pts[p1][1], pts[p3][2]-pts[p1][2]} };
                  dualmat = { {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])},
                              {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])},
                              {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                      
                  
                  break;
               }
               case 2 :
               {
                  prim_inds[0] = 4*vol+0; prim_inds[1] = 4*vol+1; prim_inds[2] = 4*vol+3;
                  dual_inds[0] = dual_index[p2]+global_node_offset_pairs[elem_edges[vol][1]][1];
                  dual_inds[1] = dual_index[p2]+global_node_offset_pairs[elem_edges[vol][2]][1];
                  dual_inds[2] = dual_index[p2]+global_node_offset_pairs[elem_edges[vol][5]][0];
                  
                  primmat = { {pts[p2][0]-pts[p0][0], pts[p2][1]-pts[p0][1], pts[p2][2]-pts[p0][2]},
                              {pts[p2][0]-pts[p1][0], pts[p2][1]-pts[p1][1], pts[p2][2]-pts[p1][2]},
                              {pts[p3][0]-pts[p2][0], pts[p3][1]-pts[p2][1], pts[p3][2]-pts[p2][2]} };
                  dualmat = { {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])},
                              {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                              {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                  break;
               }
               default :
               {
                  prim_inds[0] = 4*vol+0; prim_inds[1] = 4*vol+1; prim_inds[2] = 4*vol+2;
                  dual_inds[0] = dual_index[p3]+global_node_offset_pairs[elem_edges[vol][3]][1];
                  dual_inds[1] = dual_index[p3]+global_node_offset_pairs[elem_edges[vol][4]][1];
                  dual_inds[2] = dual_index[p3]+global_node_offset_pairs[elem_edges[vol][5]][1];
                  
                  primmat = { {pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1], pts[p3][2]-pts[p0][2]},
                              {pts[p3][0]-pts[p1][0], pts[p3][1]-pts[p1][1], pts[p3][2]-pts[p1][2]},
                              {pts[p3][0]-pts[p2][0], pts[p3][1]-pts[p2][1], pts[p3][2]-pts[p2][2]} };
                  dualmat = { {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])}, 
                              {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                              {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])} };
                  break;
               }
            }
            
            blaze::invert(primmat); //blaze::trans(primmat);
            blaze::invert(dualmat); //blaze::trans(dualmat);
            
            //~ std::cout << dualmat << std::endl;
            
            num_ele[0] = local_U[prim_inds[0]]; num_ele[1] = local_U[prim_inds[1]]; num_ele[2] = local_U[prim_inds[2]]; 
            num_mag[0] = local_F[dual_inds[0]]; num_mag[1] = local_F[dual_inds[1]]; num_mag[2] = local_F[dual_inds[2]];
            
            num_ele = blaze::eval(dualmat*num_ele);
            num_mag = blaze::eval(2*primmat*num_mag);
            
            if (abscissa.size()>0)
            {
               ordinata[0].push_back(num_ele[0]);
               ordinata[1].push_back(num_ele[1]);
               ordinata[2].push_back(num_ele[2]);
               ordinata[3].push_back(num_mag[0]);
               ordinata[4].push_back(num_mag[1]);
               ordinata[5].push_back(num_mag[2]);
            }
            else
            {
               exs[p].push_back(num_ele[0]);
               eys[p].push_back(num_ele[1]);
               ezs[p].push_back(num_ele[2]);
               hxs[p].push_back(num_mag[0]);
               hys[p].push_back(num_mag[1]);
               hzs[p].push_back(num_mag[2]);
               printed_times.push_back(t);
            }
            
            //~ auto ppt = (*op).Probe(p);
            
            //~ probe_file  << std::setw(15) << ppt[0]     << std::setw(15) << ppt[1]     << std::setw(15) << ppt[2] 
                        //~ << std::setw(15) << t
                        //~ << std::setw(15) << num_ele[0] << std::setw(15) << num_ele[1] << std::setw(15) << num_ele[2]
                        //~ << std::setw(15) << num_mag[0] << std::setw(15) << num_mag[1] << std::setw(15) << num_mag[2];
            
            
            //~ if (have_analytic)
            //~ {
               //~ auto anal_value1 = std::make_pair<Eigen::Vector3d,Eigen::Vector3d>
                                  //~ (Eigen::Vector3d({0,0,0}),Eigen::Vector3d({0,0,0}));
               //~ auto anal_value2 = std::make_pair<Eigen::Vector3d,Eigen::Vector3d>
                                  //~ (Eigen::Vector3d({0,0,0}),Eigen::Vector3d({0,0,0}));
                                  
               //~ auto analsrctype = Sources[*(Simulations[current_simulation].Src().begin())].Type();
               
               //~ if (analsrctype == "h")
               //~ {
                  //~ anal_value1 = analytic_value_excite_h(stp,
                                //~ Materials[vol_material[probe_elem[p]]].Sigma()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Epsilon()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Mu()(0,0),this->excitation_freq);
                  //~ anal_value2 = analytic_value_excite_h(stp2,Materials[vol_material[probe_elem[p]]].Sigma()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Epsilon()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Mu()(0,0),this->excitation_freq);
               //~ }
               //~ else if (analsrctype == "e")
               //~ {
                  //~ anal_value1 = analytic_value_old(stp,Materials[vol_material[probe_elem[p]]].Sigma()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Epsilon()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Mu()(0,0),this->excitation_freq);
                  //~ anal_value2 = analytic_value_old(stp2,Materials[vol_material[probe_elem[p]]].Sigma()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Epsilon()(0,0),
                                //~ Materials[vol_material[probe_elem[p]]].Mu()(0,0),this->excitation_freq);
               //~ }

               //~ analytic_e_value_vector[p].push_back(anal_value1);
               //~ analytic_h_value_vector[p].push_back(anal_value2);
            //~ }
            
            //~ probe_file  << std::endl;
         }
         
         // Clear previous plot
         plt::clf();
         // Plot a line whose name will show up as "log(x)" in the legend.
         if (abscissa.size()>0)
         {
            plt::named_plot("Ex", abscissa, ordinata[0]);
            plt::named_plot("Ey", abscissa, ordinata[1]);
            plt::named_plot("Ez", abscissa, ordinata[2]);
            plt::named_plot("Hx", abscissa, ordinata[3]);
            plt::named_plot("Hy", abscissa, ordinata[4]);
            plt::named_plot("Hz", abscissa, ordinata[5]);
            plt::xlim(abscissa[0],abscissa[probe_elem.size()-1]);
         }
         else
         {
            for (uint32_t p=0; p<probe_elem.size(); ++p)
            {
               plt::named_plot("Ex", printed_times, exs[p]);
               plt::named_plot("Ey", printed_times, eys[p]);
               plt::named_plot("Ez", printed_times, ezs[p]);
               plt::named_plot("Hx", printed_times, hxs[p]);
               plt::named_plot("Hy", printed_times, hys[p]);
               plt::named_plot("Hz", printed_times, hzs[p]);
            }
            const double xrightlim(printed_times[i] + t_step);
            plt::xlim(0.0, xrightlim);
         }
         // Set x-axis to interval [0,1000000]
         //~ plt::ylim(-1.5, 1.5);

         // Add graph title
         //~ plt::title("Sample figure");
         // Enable legend and grid
         plt::legend();
         plt::grid(true);
         // Display plot continuously
         plt::pause(0.01);
         
         //~ std::stringstream fname;
         //~ fname << (*op).Name() << i << ".png";
         //~ plt::save(fname.str());

         //~ probe_file.close();
      }
      //~ else
      //~ {
         //~ auto op = &Outputs[Simulations[current_simulation].Output()];
         //~ std::cout << "output goes in file: " << (*op).Name().c_str() << std::endl;
      //~ }
      return;
   }

   bool SameSide(Blaze3Vec v1, Blaze3Vec v2, Blaze3Vec v3, Blaze3Vec v4, Blaze3Vec p)
   {
       Blaze3Vec normal = ( (v2 - v1) % (v3 - v1) );
       double dotV4 = ( normal,(v4 - v1) );
       double dotP  = ( normal,(p - v1) );
       return ( (dotV4>-1e-6 && dotP >-1e-6) || (dotV4 < 1e-6 && dotP < 1e-6) );
       //~ return (std::signbit(dotV4) == std::signbit(dotP));
   }

   std::pair<uint32_t,uint8_t> FindProbe(const std::array<double,3UL>& pp)
   {
      Blaze3Vec pv({pp[0],pp[1],pp[2]});
      geom_tol = 1e-12;
      uint32_t v = root;
      std::vector<uint32_t> colour(elements.size()), p_queue;
      uint32_t k, qtop;
      uint32_t  curr_n;
      
      
      p_queue.push_back(root);
      colour[0]++;
      k=0;
      // std::cout << pv << std::endl;
      while ( k<p_queue.size() )
      {
         qtop=p_queue[k];
         
         v = qtop;
         
         //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         std::array<double,3UL> p1(pts[std::get<0>(elements[v])]);
         std::array<double,3UL> p2(pts[std::get<1>(elements[v])]);
         std::array<double,3UL> p3(pts[std::get<2>(elements[v])]);
         std::array<double,3UL> p4(pts[std::get<3>(elements[v])]);
         //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         

         auto v1 = Blaze3Vec({p1[0],p1[1],p1[2]});
         auto v2 = Blaze3Vec({p2[0],p2[1],p2[2]});
         auto v3 = Blaze3Vec({p3[0],p3[1],p3[2]});
         auto v4 = Blaze3Vec({p4[0],p4[1],p4[2]});

         auto check1 = SameSide(v1, v2, v3, v4, pv);
         auto check2 = SameSide(v2, v3, v4, v1, pv);
         auto check3 = SameSide(v3, v4, v1, v2, pv);
         auto check4 = SameSide(v4, v1, v2, v3, pv);
          
         if (check1 && check2 && check3 && check4)
         {
            root = v; //save root to exploit locality in next probe research
            
            double dist = blaze::norm(v1 - pv);
            uint8_t i = 0;
            if (blaze::norm(v2-pv) < dist)
            {
               i = 1;
               dist = blaze::norm(v2-pv);
            }
            if (blaze::norm(v3-pv) < dist)
            {
               i = 2;
               dist = blaze::norm(v3-pv);
            }
            if (blaze::norm(v4-pv) < dist)
               i = 3;
            
            //~ std::cout << v << std::endl;
            return std::make_pair(v,i);
         }

         for (auto curr_e : this->elem_facets[qtop])
         {   
            curr_n=*(this->facet_elems[curr_e].begin());

            if (curr_n==qtop)
               curr_n=*(this->facet_elems[curr_e].rbegin());

            if (colour[curr_n]==0)
            {
               p_queue.push_back(curr_n);
               colour[curr_n]++;
            }
         }
         colour[qtop]++;
         k++;
      }
      
      return std::make_pair(elements.size(),0);
   }
   
   bool ReadMesh(Mesh& msh)
   {
      timecounter t_preproc;
      //~ std::cout << "Loading mesh... ";
      //~ std::cout.flush();
      //~ t_preproc.tic();
      
      timecounter tc, tctot;
      double scale = msh.Scale();
      std::string input_mesh_file = msh.FileName();
      
      if (msh.IsLoaded())
         return true;
      else
         msh.Switch();
      
      /* Open file */
      if (input_mesh_file.size() == 0)
      {
         std::cout << "Invalid mesh file name" << std::endl;
         return false;
      }

      if (msh.GetMesher() == "gmsh")
      {
         if (!ConvertFromGMSH(input_mesh_file))
            return false;
         input_mesh_file = GetBaseFilename(input_mesh_file.c_str()) + std::string(".mesh");
      }
      
      uint32_t   lines, linecount;
      
      mapped_file mf(input_mesh_file);
      
      tctot.tic();
      
      /************************ Read points ************************/
      linecount = 0;
      
      xmin = 1e12;
      xmax =-1e12;
      ymin = 1e12;
      ymax =-1e12;
      zmin = 1e12;
      zmax =-1e12;
      max_edge_len = 0;
      min_height = 1e12;
      
      const char *data = mf.mem();
      
      if ( ((*data) < '0') || ((*data) > '9') )
      {
         // std::cout << "This was not supposed to happen! Invalid tetrahedral mesh input file" << std::endl;
         MyThrow(0,mesh_unknown_type);
      }
      
      char *endptr;
      
      lines = strtot<uint32_t>(data, &endptr);
      pts.reserve(lines);
      //~ std::vector<std::vector<uint32_t>> n_e;
      tc.tic();
      
      while (linecount < lines)
      {

         auto t = parser::read_point_line<double>(endptr, &endptr);
         
         
         std::array<double,3UL> point{scale*std::get<0>(t),scale*std::get<1>(t),scale*std::get<2>(t)};
         //~ point = {std::get<0>(t),std::get<1>(t),std::get<2>(t)};
         pts.push_back(point);
         n_e.push_back(std::vector<uint32_t>({}));

         if (point[0]<xmin)
            xmin = point[0];
         if (point[0]>xmax)
            xmax = point[0];
         if (point[1]<ymin)
            ymin = point[1];
         if (point[1]>ymax)
            ymax = point[1];
         if (point[2]<zmin)
            zmin = point[2];
         if (point[2]>zmax)
            zmax = point[2];
            
         ++linecount;
      }
      tc.toc();
      
      Blaze3Vec MainDiag({xmax-xmin,ymax-ymin,zmax-zmin});
      
      this->geom_tol = scale*blaze::norm(MainDiag);
      // std::cout << "Reading points: " << linecount;
      // std::cout << "/" << lines << " - " << tc << " seconds" << std::endl;
      
      /************************ Read tetrahedra ************************/
      linecount = 0;
      
      lines = strtot<uint32_t>(endptr, &endptr);
      //~ std::vector< tm_tuple > temp_tet;
      elements.reserve(lines);
      //~ elem_mat_lab.reserve(lines);
      BlazeDVec local_U(4*lines);
      BlazeDMat ETA(4UL*lines,4UL), dual_provisory(12UL*lines,3UL);
      local_U = 0;
      this->local_U = std::move(local_U);
      //~ local_U.reserve(lines);
      
      max_circum_diameter = 0;
      
      //~ std::cout << std::endl << "***************************************" << std::endl;
      //~ tc.tic();
      
      std::vector<std::pair<std::tuple<uint32_t,uint32_t>,uint32_t>> edge_list;
      std::vector<std::pair<std::tuple<uint32_t,uint32_t,uint32_t>,uint32_t>> facet_list;
      std::vector<std::vector<uint32_t>> facet_elems;
      std::vector<std::vector<std::pair<uint32_t,uint32_t>>> node_vols_w_pos(pts.size());
      
      //~ std::vector<std::vector<uint32_t>> edge_elems;
      
      
      std::vector<std::vector<Blaze3Mat>> pnt_mats;
      while (linecount < lines)
      {
         auto t = parser::read_tetrahedron_line<uint32_t>(endptr, &endptr);
         
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
                     
         elements.push_back(volume_type(p0, p1, p2, p3));
         //~ elem_mat_lab.push_back(d);
         
         std::array<double,3UL> vb0({0.25*(pts[p0][0]+pts[p1][0]+pts[p2][0]+pts[p3][0]),
                                     0.25*(pts[p0][1]+pts[p1][1]+pts[p2][1]+pts[p3][1]),
                                     0.25*(pts[p0][2]+pts[p1][2]+pts[p2][2]+pts[p3][2])});
         std::array<double,3UL> fb0({(pts[p1][0]+pts[p2][0]+pts[p3][0])/double(3),
                                     (pts[p1][1]+pts[p2][1]+pts[p3][1])/double(3),
                                     (pts[p1][2]+pts[p2][2]+pts[p3][2])/double(3)});
         std::array<double,3UL> fb1({(pts[p0][0]+pts[p2][0]+pts[p3][0])/double(3),
                                     (pts[p0][1]+pts[p2][1]+pts[p3][1])/double(3),
                                     (pts[p0][2]+pts[p2][2]+pts[p3][2])/double(3)});
         std::array<double,3UL> fb2({(pts[p0][0]+pts[p1][0]+pts[p3][0])/double(3),
                                     (pts[p0][1]+pts[p1][1]+pts[p3][1])/double(3),
                                     (pts[p0][2]+pts[p1][2]+pts[p3][2])/double(3)});
         std::array<double,3UL> fb3({(pts[p0][0]+pts[p1][0]+pts[p2][0])/double(3),
                                     (pts[p0][1]+pts[p1][1]+pts[p2][1])/double(3),
                                     (pts[p0][2]+pts[p1][2]+pts[p2][2])/double(3)});
         
         node_vols_w_pos[p0].push_back(std::make_pair(linecount,0));
         node_vols_w_pos[p1].push_back(std::make_pair(linecount,1));
         node_vols_w_pos[p2].push_back(std::make_pair(linecount,2));
         node_vols_w_pos[p3].push_back(std::make_pair(linecount,3));
         
         
         //~ std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         
         //~ std::cout << vb0;
         //~ std::cout << fb0 << fb1 << fb2 << fb3;
         
         //~ Blaze3Mat prim_0, prim_1, prim_2, prim_3;
         auto prim_0 = blaze::submatrix<blaze::aligned>(dual_provisory,12UL*linecount+0UL,0UL,3UL,3UL);
         auto prim_1 = blaze::submatrix<blaze::aligned>(dual_provisory,12UL*linecount+3UL,0UL,3UL,3UL);
         auto prim_2 = blaze::submatrix<blaze::aligned>(dual_provisory,12UL*linecount+6UL,0UL,3UL,3UL);
         auto prim_3 = blaze::submatrix<blaze::aligned>(dual_provisory,12UL*linecount+9UL,0UL,3UL,3UL);
         
         Blaze3Mat dual_0, dual_1, dual_2, dual_3;
         
         prim_0 = { {pts[p1][0]-pts[p0][0], pts[p1][1]-pts[p0][1], pts[p1][2]-pts[p0][2]},
                    {pts[p2][0]-pts[p0][0], pts[p2][1]-pts[p0][1], pts[p2][2]-pts[p0][2]},
                    {pts[p3][0]-pts[p0][0], pts[p3][1]-pts[p0][1], pts[p3][2]-pts[p0][2]} };
         
         prim_1 = { {-pts[p1][0]+pts[p0][0], -pts[p1][1]+pts[p0][1], -pts[p1][2]+pts[p0][2]},
                    {pts[p2][0]-pts[p1][0], pts[p2][1]-pts[p1][1], pts[p2][2]-pts[p1][2]},
                    {pts[p3][0]-pts[p1][0], pts[p3][1]-pts[p1][1], pts[p3][2]-pts[p1][2]} };
         
         prim_2 = { {-pts[p2][0]+pts[p0][0], -pts[p2][1]+pts[p0][1], -pts[p2][2]+pts[p0][2]},
                    {-pts[p2][0]+pts[p1][0], -pts[p2][1]+pts[p1][1], -pts[p2][2]+pts[p1][2]},
                    {pts[p3][0]-pts[p2][0], pts[p3][1]-pts[p2][1], pts[p3][2]-pts[p2][2]} };
                    
         prim_3 = { {-pts[p3][0]+pts[p0][0], -pts[p3][1]+pts[p0][1], -pts[p3][2]+pts[p0][2]},
                    {-pts[p3][0]+pts[p1][0], -pts[p3][1]+pts[p1][1], -pts[p3][2]+pts[p1][2]},
                    {-pts[p3][0]+pts[p2][0], -pts[p3][1]+pts[p2][1], -pts[p3][2]+pts[p2][2]} };


         auto volume = det(prim_0);
         std::array<double,4> coeffs({1,-1,1,-1});
         if (volume<0)
         {
            coeffs[0]=-coeffs[0];
            coeffs[1]=-coeffs[1];
            coeffs[2]=-coeffs[2];
            coeffs[3]=-coeffs[3];
         }
         
         //~ std::cout << blaze::row<1UL>(prim_1) << std::endl;
         //~ std::cout << blaze::row<2UL>(prim_1) << std::endl;
         blaze::StaticVector<double,3UL,blaze::rowVector> trig_0 = (blaze::row<1UL>(prim_1) % blaze::row<2UL>(prim_1));
         blaze::StaticVector<double,3UL,blaze::rowVector> trig_1 = (blaze::row<0UL>(prim_2) % blaze::row<0UL>(prim_3));
         blaze::StaticVector<double,3UL,blaze::rowVector> trig_2 = (blaze::row<0UL>(prim_0) % blaze::row<0UL>(prim_3));
         blaze::StaticVector<double,3UL,blaze::rowVector> trig_3 = (blaze::row<0UL>(prim_0) % blaze::row<1UL>(prim_0));
         
         std::cout << std::endl;
         
         const double norm_0 = blaze::norm(trig_0);
         const double norm_1 = blaze::norm(trig_1);
         const double norm_2 = blaze::norm(trig_2);
         const double norm_3 = blaze::norm(trig_3);
         
         double height_0 = std::fabs(volume)/std::fabs(norm_0);
         double height_1 = std::fabs(volume)/std::fabs(norm_1);
         double height_2 = std::fabs(volume)/std::fabs(norm_2);
         double height_3 = std::fabs(volume)/std::fabs(norm_3);
         
         if (height_0 < min_height)
            min_height = height_0;
         if (height_1 < min_height)
            min_height = height_1;
         if (height_2 < min_height)
            min_height = height_2;
         if (height_3 < min_height)
            min_height = height_3;
         
         volume = std::fabs(volume)/double(6); // get the actual tet's volume
            
         dual_0 = {  {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                     {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])},
                     {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                      
         dual_1 = {  {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])},
                     {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])},
                     {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                      
         dual_2 = {  {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])},
                     {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                     {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                      
         dual_3 = {  {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])}, 
                     {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                     {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])} };

         //~ std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         
         //~ std::cout << "dual matrix test:" << std::endl << dual_0;
         blaze::invert(prim_0); //blaze::trans(prim_0);
         blaze::invert(prim_1); //blaze::trans(prim_1);
         blaze::invert(prim_2); //blaze::trans(prim_2);
         blaze::invert(prim_3); //blaze::trans(prim_3);
         blaze::invert(dual_0);
         blaze::invert(dual_1);
         blaze::invert(dual_2);
         blaze::invert(dual_3);
         
         prim_0 = volume*blaze::eval(blaze::trans(prim_0)*blaze::eval(Materials[d].Mu()*prim_0));
         prim_1 = volume*blaze::eval(blaze::trans(prim_1)*blaze::eval(Materials[d].Mu()*prim_1));
         prim_2 = volume*blaze::eval(blaze::trans(prim_2)*blaze::eval(Materials[d].Mu()*prim_2));
         prim_3 = volume*blaze::eval(blaze::trans(prim_3)*blaze::eval(Materials[d].Mu()*prim_3));
         
         dual_0 = 0.25*volume*blaze::eval(blaze::trans(dual_0)*blaze::eval(Materials[d].Epsilon()*dual_0));
         dual_1 = 0.25*volume*blaze::eval(blaze::trans(dual_1)*blaze::eval(Materials[d].Epsilon()*dual_1));
         dual_2 = 0.25*volume*blaze::eval(blaze::trans(dual_2)*blaze::eval(Materials[d].Epsilon()*dual_2));
         dual_3 = 0.25*volume*blaze::eval(blaze::trans(dual_3)*blaze::eval(Materials[d].Epsilon()*dual_3));
         
         
         auto prim_local_MM = blaze::submatrix(ETA,4*linecount,0UL,4UL,4UL);
         prim_local_MM = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
         
         prim_local_MM(1,1) += dual_0(0,0); prim_local_MM(1,2) += dual_0(0,1); prim_local_MM(1,3) += dual_0(0,2);
         prim_local_MM(2,1) += dual_0(1,0); prim_local_MM(2,2) += dual_0(1,1); prim_local_MM(2,3) += dual_0(1,2);
         prim_local_MM(3,1) += dual_0(2,0); prim_local_MM(3,2) += dual_0(2,1); prim_local_MM(3,3) += dual_0(2,2);
         
         prim_local_MM(0,0) += dual_1(0,0); prim_local_MM(0,2) += dual_1(0,1); prim_local_MM(0,3) += dual_1(0,2);
         prim_local_MM(2,0) += dual_1(1,0); prim_local_MM(2,2) += dual_1(1,1); prim_local_MM(2,3) += dual_1(1,2);
         prim_local_MM(3,0) += dual_1(2,0); prim_local_MM(3,2) += dual_1(2,1); prim_local_MM(3,3) += dual_1(2,2);
         
         prim_local_MM(0,0) += dual_2(0,0); prim_local_MM(0,1) += dual_2(0,1); prim_local_MM(0,3) += dual_2(0,2);
         prim_local_MM(1,0) += dual_2(1,0); prim_local_MM(1,1) += dual_2(1,1); prim_local_MM(1,3) += dual_2(1,2);
         prim_local_MM(3,0) += dual_2(2,0); prim_local_MM(3,1) += dual_2(2,1); prim_local_MM(3,3) += dual_2(2,2);
         
         prim_local_MM(0,0) += dual_3(0,0); prim_local_MM(0,1) += dual_3(0,1); prim_local_MM(0,2) += dual_3(0,2);
         prim_local_MM(1,0) += dual_3(1,0); prim_local_MM(1,1) += dual_3(1,1); prim_local_MM(1,2) += dual_3(1,2);
         prim_local_MM(2,0) += dual_3(2,0); prim_local_MM(2,1) += dual_3(2,1); prim_local_MM(2,2) += dual_3(2,2);
         
         blaze::StaticMatrix<double,4UL,3UL> dual_edges, primal_faces;
         
         dual_edges = {  {coeffs[0]*(fb0[0]-vb0[0]),coeffs[0]*(fb0[1]-vb0[1]),coeffs[0]*(fb0[2]-vb0[2])},
                         {coeffs[1]*(fb1[0]-vb0[0]),coeffs[1]*(fb1[1]-vb0[1]),coeffs[1]*(fb1[2]-vb0[2])},
                         {coeffs[2]*(fb2[0]-vb0[0]),coeffs[2]*(fb2[1]-vb0[1]),coeffs[2]*(fb2[2]-vb0[2])},
                         {coeffs[3]*(fb3[0]-vb0[0]),coeffs[3]*(fb3[1]-vb0[1]),coeffs[3]*(fb3[2]-vb0[2])} };
                         
         primal_faces = {  {0.5*trig_0[0],0.5*trig_0[1],0.5*trig_0[2]},
                           {0.5*trig_1[0],0.5*trig_1[1],0.5*trig_1[2]},
                           {0.5*trig_2[0],0.5*trig_2[1],0.5*trig_2[2]},
                           {0.5*trig_3[0],0.5*trig_3[1],0.5*trig_3[2]}  };
   
         
         //~ auto pippo = prim_local_MM*dual_edges;
         
         //~ std::cout << "consistency check:" << std::endl;
         //~ std::cout << pippo << std::endl;
         //~ std::cout << primal_faces << std::endl;

         blaze::invert<blaze::asSymmetric>(prim_local_MM);
         
         std::array<std::tuple<uint32_t,uint32_t,uint32_t>,4> local_facets;
         local_facets[0] = std::make_tuple(p0,p1,p2);
         local_facets[1] = std::make_tuple(p0,p1,p3);
         local_facets[2] = std::make_tuple(p0,p2,p3);
         local_facets[3] = std::make_tuple(p1,p2,p3);

         std::array<std::tuple<uint32_t,uint32_t>,6> local_edges;
         local_edges[0] = std::make_tuple(p0,p1);
         local_edges[1] = std::make_tuple(p0,p2);
         local_edges[2] = std::make_tuple(p1,p2);
         local_edges[3] = std::make_tuple(p0,p3);
         local_edges[4] = std::make_tuple(p1,p3);
         local_edges[5] = std::make_tuple(p2,p3);
         
         for (auto lf : local_facets)
            facet_list.push_back(std::make_pair(lf,linecount));
         
         for (auto le : local_edges)
            edge_list.push_back(std::make_pair(le,linecount));
         
         linecount++;
      }
      //~ std::cout << "---------------------------------------------------------------" << std::endl;
      
      std::vector<std::vector<uint32_t>> elem_edges(lines);
      std::vector<std::vector<uint32_t>> elem_facets(lines);
      
      /************************ Read boundary surfaces ************************/
      linecount = 0;
      //~ auto num_of_tets=lines;
      lines = strtot<uint32_t>(endptr, &endptr);
      std::vector<std::set<uint32_t>> bndnd(pts.size());
      
      while (linecount < lines)
      {
         
         // if ( (linecount%50000) == 0 )
         // {
            // std::cout << "Reading triangles: " << linecount;
            // std::cout << "/" << lines << "\r";
            // std::cout.flush();
         // }
         
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
         //~ surface_type   tri( p0, p1, p2 );
         
         bndnd[p0].insert(bid);
         bndnd[p1].insert(bid);
         bndnd[p2].insert(bid);
         
         
         ++linecount;
      }
      
      mf.close();
      this->bndnd = std::move(bndnd);
      
      /************************ Generate Data Structures ************************/
      struct {
         bool operator()(const std::pair<std::tuple<uint32_t,uint32_t>,uint32_t> t1,
                         const std::pair<std::tuple<uint32_t,uint32_t>,uint32_t> t2)
         {
            return (t1.first < t2.first);
         }
      } compedges;
      
      struct {
         bool operator()(const std::pair<std::tuple<uint32_t,uint32_t,uint32_t>,uint32_t> t1,
                         const std::pair<std::tuple<uint32_t,uint32_t,uint32_t>,uint32_t> t2)
         {
            return (t1.first < t2.first);
         }
      } compfacets;
      

      //~ std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      
      //~ std::mutex max_mutex;
      
      std::thread thread1( [&] 
      { 
         std::sort(edge_list.begin(),edge_list.end(),compedges);
         auto ptr1 = edge_list.begin();
         auto ptr2 = edge_list.begin();
         ptr2++;
         
         uint32_t j=0;
         edges.push_back((*ptr1).first);
         
         //~ edge_elems.push_back(std::vector<uint32_t>({}));
         //~ edge_elems[j].push_back((*ptr1).second);
         
         auto n0 = std::get<0>((*ptr1).first);
         auto n1 = std::get<1>((*ptr1).first);
         global_node_offset_pairs.push_back(std::array<uint32_t,2>({n_e[n0].size(),n_e[n1].size()}));
         
         n_e[n0].push_back(j);
         n_e[n1].push_back(j);
         elem_edges[(*ptr1).second].push_back(j);
         
         while ( ptr2 != edge_list.end())
         {
            auto vv = (*ptr2).second;
            if ((*ptr2).first != (*ptr1).first)
            {
               edges.push_back((*ptr2).first);
               
               ++j;
               
               //~ edge_elems.push_back(std::vector<uint32_t>({}));
               //~ edge_elems[++j].push_back(vv);
               
               auto n0 = std::get<0>((*ptr2).first);
               auto n1 = std::get<1>((*ptr2).first);
               
               global_node_offset_pairs.push_back(std::array<uint32_t,2>({n_e[n0].size(),n_e[n1].size()}));
               
               n_e[n0].push_back(j);
               n_e[n1].push_back(j);
               
            }
            
            elem_edges[vv].push_back(j);
            if (elem_edges[vv].size()==6)
               std::swap(elem_edges[vv][2],elem_edges[vv][3]);
            
            ++ptr1;
            ++ptr2;
         }
      });
      
      std::thread thread2( [&] 
      { 
         std::sort(facet_list.begin(),facet_list.end(),compfacets);
         auto ptr1 = facet_list.begin();
         auto ptr2 = facet_list.begin();
         ptr2++;
         
         uint32_t j=0;
         
         facets.push_back((*ptr1).first);
         facet_elems.push_back(std::vector<uint32_t>({}));
         facet_elems[j].push_back((*ptr1).second);
         elem_facets[(*ptr1).second].push_back(j);
         while ( ptr2 != facet_list.end())
         {
            if ((*ptr2).first != (*ptr1).first)
            {
               facets.push_back((*ptr2).first);
               facet_elems.push_back(std::vector<uint32_t>({}));
               ++j;
            }
            
            facet_elems[j].push_back((*ptr2).second);
            elem_facets[(*ptr2).second].push_back(j);
            
            if (elem_facets[(*ptr2).second].size()==4)
            {
               auto vv = (*ptr2).second;
               std::swap(elem_facets[vv][0],elem_facets[vv][3]);
               std::swap(elem_facets[vv][1],elem_facets[vv][2]);
            }
            
            ++ptr1;
            ++ptr2;
         }
      });
      
      thread1.join();
      thread2.join();
      
      uint32_t tot_local_dual=0;
      
      for (uint32_t j=0;j<pts.size();++j)
      {
         blaze::HybridMatrix<double,25UL,25UL> local_M(16,16); //magic number in node-edge incidence
         local_M = 0;
         
         for (auto vol : node_vols_w_pos[j])
         {
            //~ auto sprim = pnt_mats[vol.first][vol.second];
            auto sprim = blaze::submatrix(dual_provisory,12UL*vol.first+3UL*vol.second,0UL,3UL,3UL);
            uint32_t i1,i2,i3;
            
            switch (vol.second)
            {
               case 0 :
               {
                  i1 = global_node_offset_pairs[elem_edges[vol.first][0]][0];//.second;
                  i2 = global_node_offset_pairs[elem_edges[vol.first][1]][0];//.second;
                  i3 = global_node_offset_pairs[elem_edges[vol.first][3]][0];//.second;
                  break;
               }
               case 1 :
               {
                  i1 = global_node_offset_pairs[elem_edges[vol.first][0]][1];//.second;
                  i2 = global_node_offset_pairs[elem_edges[vol.first][2]][0];//.second;
                  i3 = global_node_offset_pairs[elem_edges[vol.first][4]][0];//.second;
                  break;
               }
               case 2 :
               {
                  i1 = global_node_offset_pairs[elem_edges[vol.first][1]][1];//.second;
                  i2 = global_node_offset_pairs[elem_edges[vol.first][2]][1];//.second;
                  i3 = global_node_offset_pairs[elem_edges[vol.first][5]][0];//.second;
                  break;
               }
               default:
               {
                  i1 = global_node_offset_pairs[elem_edges[vol.first][3]][1];//.second;
                  i2 = global_node_offset_pairs[elem_edges[vol.first][4]][1];//.second;
                  i3 = global_node_offset_pairs[elem_edges[vol.first][5]][1];//.second;
                  break;
               }
            }
            
            local_M(i1,i1) += sprim(0,0); local_M(i1,i2) += sprim(0,1); local_M(i1,i2) += sprim(0,2);
            local_M(i2,i1) += sprim(1,0); local_M(i2,i2) += sprim(1,1); local_M(i2,i3) += sprim(1,2);
            local_M(i3,i1) += sprim(2,0); local_M(i3,i2) += sprim(2,1); local_M(i3,i3) += sprim(2,2);
         }
         
         //~ for (uint8_t k=0;k<n_e[j].size();++k)
            //~ for (uint8_t l=k+1;l<n_e[j].size();++l)
               //~ local_M(k,l) = local_M(l,k);
               
         
         //padding
         if (n_e[j].size()<16)
            for(uint8_t k=n_e[j].size();k<16;++k)
               local_M(k,k) = 1.0;
         
         //~ std::cout << local_M << std::endl;
         
         //~ std::cout << __FILE__ << ":" << __LINE__ << std::endl;
         
         blaze::invert<blaze::asSymmetric>(local_M);
         
         //~ std::cout << "Full: " << std::endl;
         //~ std::cout << local_M << std::endl;
         
         BlazeSymSpMat dual_M(n_e[j].size());
         dual_M.reserve(n_e[j].size()*n_e[j].size());
         
         for (uint32_t k=0;k<n_e[j].size();++k)
            dual_M.reserve(k,n_e[j].size());
         
         for (uint32_t k=0;k<n_e[j].size();++k)
         {
            dual_M.append(k,k,local_M(k,k));
            
            for (uint32_t l=k+1;l<n_e[j].size();++l)
               if (std::fabs(local_M(l,k))> 1e-14)
                  dual_M.append(k,l,local_M(k,l));
            
            dual_M.finalize( k );
         }
         
         //~ std::cout << "Sparse: " << std::endl;
         //~ std::cout << dual_M << std::endl;
         this->dual_index.push_back(tot_local_dual);
         tot_local_dual += n_e[j].size();
        
         //~ local_F.push_back(lfl);
         this->NU.push_back(dual_M);
      }

      
      
      BlazeDVec local_F(tot_local_dual);
      local_F = 0;
      this->local_F = std::move(local_F);
      this->ETA = std::move(ETA);
      
      this->elem_facets = std::move(elem_facets);
      this->elem_edges  = std::move(elem_edges);
      this->facet_elems = std::move(facet_elems);
      //~ this->edge_elems  = std::move(edge_elems);
      
      return true;
   }

   std::map<uint32_t,Output>                                Outputs;
   std::map<uint32_t,Simulation>                            Simulations;
   std::map<uint32_t,Source>                                Sources;
   std::map<uint32_t,Material>                              Materials;
   std::map<uint32_t,BoundaryCondition>                     BCs;
   std::map<uint32_t,Mesh>                                  Meshes;
   std::map<uint32_t,Solid>                                 Solids;
   std::map<uint32_t,Refinement>                            Refinements;
   std::vector<uint32_t>                                    BndToBC;
   uint32_t                                                 method_line;
   
   private:
   uint32_t    input_line, h_mat_fill_in, n_mat_fill_in, i_silo;
   uint32_t    loaded_mesh_label, current_simulation, root;
   bool        possibly_unstable;
   double      t_step, min_h, average_diameter, max_circum_diameter, excitation_freq;
   double      xmin,xmax,ymin,ymax,zmin,zmax,min_height,max_edge_len,geom_tol;
   Duration    simulation_time;
   
   std::vector<std::array<double,3UL>>                            pts;
   BlazeDVec                                                      local_U, local_F;
   std::vector<BlazeSymSpMat>                                     NU;
   BlazeDMat                                                      ETA;
   std::vector<std::array<uint32_t,2>>                            global_node_offset_pairs;
   std::vector<std::vector<uint32_t>>                             n_e,/*edge_elems, */facet_elems,elem_edges,elem_facets;
   std::vector<std::tuple<uint32_t,uint32_t,uint32_t>>            facets;
   std::vector<std::tuple<uint32_t,uint32_t>>                     edges;
   std::vector<std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>>   elements;
   std::vector<std::pair<uint32_t,uint8_t>>                       probe_elem;
   std::vector<uint32_t>                                          /*elem_mat_lab, */dual_index,probe_index;
   std::vector<double>                                            abscissa, printed_times;
   std::vector<std::set<uint32_t>>                                bndnd;
   std::vector<std::vector<double>>                               eys,exs,ezs,hxs,hys,hzs;
   
};
