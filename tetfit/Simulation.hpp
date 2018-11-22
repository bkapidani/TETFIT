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
      method = "dga";
      solver = "none";
      courant = 1;
      hardcoded_time_step = 0;
      tol = 1e-8;
      debug_mat = false;
      analytic = "false";
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
      else if (param == "courant")
         courant = std::stod(value);
      else if (param == "step")
         hardcoded_time_step = std::stod(value);
      else if (param == "output")
         output = std::stod(value);
      else if (param == "method")
         SetMethod(input_line,value);
      else if (param == "solver")
         SetSolver(input_line,value);
      else if (param == "tolerance")
         tol = std::stod(value);
      else if (param == "analytic")
      {
         // if (value == "true")
            // analytic = true;
         // else if (value == "false")
            // analytic = false;
         // else
            // MyThrow(input_line,sim_unknown_parameter);
         analytic = value;
      }
      else if (param == "debug")
      {
         if (value == "true")
            debug_mat = true;
         else if (value == "false")
            debug_mat = false;
         else
            MyThrow(input_line,sim_unknown_parameter);

      }
      else
         MyThrow(input_line,sim_unknown_parameter);
   }
   
   const std::string& HaveAnalytic(void) const { return analytic; }
   const bool& DebugMatrices(void) const { return debug_mat; }
   const Duration& Time(void) const { return d; } 
   const uint32_t& MeshLabel(void) const { return mesh_label; }
   const uint32_t& Output(void) const { return output; }
   const std::vector<uint32_t>& Src(void) { return sources; }
   const SimMethod& Method(void) const { return method; }
   void ForceMethod(const SimMethod& m){ SetMethod(0,m); }
   const Solver& GetSolver(void) const { return solver; }
   const double& Tolerance(void) const { return tol; }
   const double& Courant(void) const { return courant; }
   const double& HardCoded_TS(void) const { return hardcoded_time_step; }
   // const Blaze3Vec& Probe(void) const { return probepoint; }
   
   private:
   Duration d;
   SimMethod method;
   Solver solver;
   bool debug_mat;
   std::string analytic;
   double tol, courant, hardcoded_time_step;
   // Blaze3Vec probepoint;
   std::vector<uint32_t> sources; //can combine multiple sources
   uint32_t mesh_label, output;
   void SetMethod(uint32_t il, SimMethod m)
   {
      if (std::find(simmethods.begin(),simmethods.end(),m) == simmethods.end())
                        MyThrow(il,sim_unknown_method);
      else 
         method = m;
   }
   
   void SetSolver(uint32_t il, Solver m)
   {
      if (std::find(solvers.begin(),solvers.end(),m) == solvers.end())
                        MyThrow(il,sim_unknown_solver);
      else 
         solver = m;
   }   
   
};
