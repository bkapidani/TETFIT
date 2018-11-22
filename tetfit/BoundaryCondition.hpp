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
 
 //file BoundaryCondition.hpp
 #include "Utilities.hpp"
 
class BoundaryCondition
{
   public:
   BoundaryCondition()
   : type("none")
   {
      thickness=0; //just for PML, otherwise unused
      dir = "none";
      is_set = false;
   }
   
   bool Set(void) { return is_set; };
   
   void SetParam(uint32_t input_line, std::string param, std::string value)
   {
      is_set = true;
      if (param == "type")
      {
         // Temporary error messages
         if (value == "pml")
         {
            // MyThrow(input_line,pml_missing);
            type = value;
            // this->val=0;
         }
         else if (value == "pmc")
         {
            type = value;
            this->val=0;
         }
         else if (value == "pec")
         {
            type = value;
            this->val=0;
         }
         else
            MyThrow(input_line,bc_unknown_type);
      }
      else if (param == "direction")
      {
         if (std::find(bc_directions.begin(),bc_directions.end(),value) == bc_directions.end())
            MyThrow(input_line,bc_unknown_direction);
         else
            this->dir = value;
      }
      else if (param == "width")
         thickness = std::stod(value);
      else if (param == "surface")
         surfs.push_back(std::stoi(value));
      else if (param == "surfaces")
      {
         auto  i = value.begin();
         if (*i != '{')
            MyThrow(input_line,coordinates_syntax);
         else
         {
            i++;
            
            while (*i != '}' && i != value.end())
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
                  surfs.push_back(std::stoi(coord));
                  
                  if (*i == ',')
                     i++;
               }
            }
            
            if (i == value.end())
               MyThrow(input_line,unbalanced_bracket);
         }
      }
      else   
         MyThrow(input_line,bc_unknown_parameter);
   }
   
   const std::string&       Type(void) { return type; }
   const double&          Width(void) { return thickness; }
   const double&         GetValue(void) { return val; }
   const std::vector<uint32_t> Surfaces(void) { return surfs; }
   const Direction&      GetDirection(void) { return dir; }
   
   private:
   std::string type;
   Direction dir;
   bool is_set;
   double thickness, val;
   std::vector<uint32_t> surfs;
};
