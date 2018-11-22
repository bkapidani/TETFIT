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
      // probepoints = Blaze3Vec({0,0,0});
      mode = "none";
      name = "simulation";
      output_period = 0;
      nextprint = 1e-15;
      delay     = 0;
      grid = false;
      error_norm = false;
      xstart = ystart = zstart = 0;
      xstop  = ystop  = zstop  = 0;
      xstep  = ystep  = zstep  = 0;
      ref_sys = 0; //cartesian coordinates
      // silo_instant = 0;
   }
   
   void Initialize(void)
   {
      nextprint = delay;
   }
   
   void SetParam(uint32_t input_line, std::string param, std::string value)
   {
      if (param == "mode")
      {
         if (std::find(outputmodes.begin(),outputmodes.end(),value) == outputmodes.end())
            MyThrow(input_line,sim_unknown_output);
         mode = value;
         //~ if (mode == "probefield")
            //~ probepoints.reserve(1000);
      }
      else if (param == "radiator")
         radiating_vol_bnd.push_back(std::stof(value));
      else if (param == "delay")
      {
         nextprint = delay = std::stod(value);
      }
      else if (param == "coordinates")
      {
         if (value == "cartesian")
            ref_sys = 0;
         else if (value == "cylindrical")
            ref_sys = 1;
         else if (value == "spherical")
            ref_sys = 2;
         else
            MyThrow(input_line,sim_invalid_coord_system);
      }
      else if (param == "probe")
      {
         std::array<double,3UL> new_probepoint;
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
                     new_probepoint[k]= std::stof(coord);
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
         // std::cout << new_probepoint << std::endl;
         probepoints.push_back(new_probepoint);
      }
      else if (param == "period")
      {
         output_period = std::stof(value);
      }
      else if (param == "xgrid")
      {
         // double xstart,xstep,xstop;
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
                     xstart = std::stof(coord);
                     k++;
                  }
                  else if (k == 1)
                  {
                     xstep = std::stof(coord);
                     k++;
                  }
                  else if (k == 2)
                  {
                     xstop = std::stof(coord);
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
         
         assert(xstart <= xstop);
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
                     ystart = std::stof(coord);
                     k++;
                  }
                  else if (k == 1)
                  {
                     ystep = std::stof(coord);
                     k++;
                  }
                  else if (k == 2)
                  {
                     ystop = std::stof(coord);
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
         
         assert(ystart <= ystop);
         assert(ystep > 0);
      }
      else if (param == "zgrid")
      {
         // double xstart,xstep,xstop;
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
                     zstart = std::stof(coord);
                     k++;
                  }
                  else if (k == 1)
                  {
                     zstep = std::stof(coord);
                     k++;
                  }
                  else if (k == 2)
                  {
                     zstop = std::stof(coord);
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
         
         assert(zstart <= zstop);
         assert(zstep > 0);
      }
      else if (param == "grid")
      {
         if (value == "on")
         {
            double xiter, yiter, ziter;
            for (ziter=zstart; ziter<=zstop; ziter+=zstep)
               for (yiter=ystart; yiter<=ystop; yiter+=ystep)
                  for (xiter=xstart; xiter<=xstop; xiter+=xstep)
                     probepoints.push_back(std::array<double,3UL>({xiter,yiter,ziter}));
         }
         else if (value == "off")
            ;
         else
            MyThrow(input_line,out_unknown_parameter);
      }
      else if (param == "error")
      {
         if (value == "on")
            error_norm = true;
         else if (value == "off")
            error_norm = false;
         else
            MyThrow(input_line,out_unknown_parameter);
      }
      else if (param == "name")
         name = value;
      else
         MyThrow(input_line,out_unknown_parameter);
   }
   
   bool AllowPrint(double t)
   {
      
      if (t >= nextprint)
      {
         nextprint+=output_period;
         return true;
      }
      else
         return false;
   }
   
   const uint32_t Nprobes(void) const { return probepoints.size(); }
   const uint8_t ReferenceFrame(void) const { return ref_sys; }
   const double& Period(void) const { return output_period; }
   // const double& Instant(void) const { return silo_instant; }
   const std::string& Name(void) const { return name; }
   const OutputMode& Mode(void) const { return mode; }
   const std::vector<uint32_t> GetRadiators(void) const {return radiating_vol_bnd; };
   const std::array<double,3UL>& Probe(uint32_t i) const
   { 
      assert(i<probepoints.size());
      return probepoints[i];
   }
   
   private:
   bool grid, error_norm;
   std::vector<uint32_t> radiating_vol_bnd;
   std::string name;
   double xstart, xstep, xstop, ystart, ystep, ystop, zstart, zstep, zstop;
   double output_period, nextprint, delay;
   std::vector<std::array<double,3UL>> probepoints;
   OutputMode mode;
   uint8_t ref_sys;
};
