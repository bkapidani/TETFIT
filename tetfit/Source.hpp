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
 
 //file Source.hpp
 #include "Utilities.hpp"
 
 class Source
{
   public:
   Source(void)
   {
      //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      is_set     = false;
      carrier    = "sin";
      prof       = "dc";
      dir         = "x";
      amp         = 0;
      freq        = 0;
      radius      = 1;
      decay      = 0;
      width       = 1;
      //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      delay       = 0;
      impedance = 50; //Ohm (only real allowed for now)
      kvec[0] = kvec[1] = kvec[2] = 0;
      location[0] = location[1] = location[2] = 0;
      //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
      st      = " ";
      bfuncs[0] = bfuncs[1] = bfuncs[2] = "cos";
      center_coords[0] = center_coords[1] = center_coords[2] = 0;
      surface_label.resize(0);
      //std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   }
   
   Source(const Source& ) = default;
   
   bool Set(void) { return is_set; };
   
   void SetParam(uint32_t input_line, std::string param, std::string value)
   {
      is_set = true;
      if (param == "direction")
      {
         if (std::find(directions.begin(),directions.end(),value) == directions.end())
            MyThrow(input_line,src_unknown_direction);
         else
            this->dir = value;
      }
      else if (param == "type")
      {
         if (std::find(sourcetypes.begin(),sourcetypes.end(),value) == sourcetypes.end())
            MyThrow(input_line,src_unknown_type);
         else
            this->st = value;
      }
      else if (param == "profile")
      {
         if (std::find(profiles.begin(),profiles.end(),value) == profiles.end())
            MyThrow(input_line,src_unknown_profile);
         else
            this->prof = value;
      }
      else if (param == "carrier")
      {
         if (std::find(carriers.begin(),carriers.end(),value) == carriers.end())
            MyThrow(input_line,src_unknown_carrier);
         else
            this->carrier = value;
      }
      else if (param == "center")
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
                     center_coords[k]= std::stod(coord);
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
      else if (param == "wavevector")
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
                  if (k < 3)
                  {
                     kvec[k]= std::stod(coord);
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
      else if (param == "location")
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
                  if (k < 3)
                  {
                     location[k]= std::stod(coord);
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
      else if (param == "mode")
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
                  if (*i != ' ')
                     coord.push_back(*i);
                  i++;
               }
               if (i == value.end())
                  MyThrow(input_line,unbalanced_bracket);
               else 
               {
                  if (k < 3)
                  {
                     
                     // if (std::find(modes.begin(),modes.end(),coord) == modes.end())
                        // MyThrow(input_line,src_unknown_bf);
                     
                     bfuncs[k] = coord;
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
      else if (param == "Xmode")
      {
         if (std::find(modes.begin(),modes.end(),value) == modes.end())
            MyThrow(input_line,src_unknown_bf);
         else
            this->bfuncs[0] = value;
      }
      else if (param == "Ymode")
      {
         if (std::find(modes.begin(),modes.end(),value) == modes.end())
            MyThrow(input_line,src_unknown_bf);
         else
            this->bfuncs[1] = value;
      }
      else if (param == "Zmode")
      {
         if (std::find(modes.begin(),modes.end(),value) == modes.end())
            MyThrow(input_line,src_unknown_bf);
         else
            this->bfuncs[2] = value;
      }
      else if (param == "amplitude")
         this->amp = std::stod(value);
      else if (param == "delay")
         this->delay = std::stod(value);
      else if (param == "impedance")
         this->impedance = std::stod(value);
      else if (param == "frequency")
         this->freq = std::stod(value);
      else if (param == "radius")
         this->radius = std::stod(value);
      else if (param == "decay")
         this->decay = std::stod(value);
      else if (param == "width")
         this->width = std::stod(value);
      else if (param == "kx")
         this->kvec[0] = std::stod(value);
      else if (param == "ky")
         this->kvec[1] = std::stod(value);
      else if (param == "kz")
         this->kvec[2] = std::stod(value);
      else if (param == "surface")
         this->surface_label.push_back(std::stoi(value));
      else
         MyThrow(input_line,src_unknown_parameter);
   }
   
   const Sourcetype&    Type(void) { return st; }
   const std::vector<uint32_t>&     Surface() { return surface_label; } 
   const Direction&    GetDirection(void) { return dir; }
   const WaveNumber&   Getkx(void) { return kvec[0]; }
   const WaveNumber&   Getky(void) { return kvec[1]; }
   const WaveNumber&   Getkz(void) { return kvec[2]; }
   const Frequency&   GetFreq(void) { return freq; }
   const Amplitude&   GetAmp(void) { return amp; }
   const Blaze3Vec& Location(void) { return location; }
   
   Blaze3Vec Compute(SpaceTimePoint p) //returns a vector, so one can use superposition of sources
   {
      double ret=amp; // if the source is dc, we are already done!
      if (st == "j")
         ret *= impedance;
      /*Time profile of source*/
      if (carrier == "sin")
         ret *= sin(2*PI*freq*(p[3]-delay)); 
      else if (carrier == "cos")
         ret *= amp*cos(2*PI*freq*(p[3]-delay));
      else if (carrier == "gaussian")
         ret *= std::exp(-std::pow(p[3]-delay,2)/std::pow(width,2));
      else if (carrier == "ricker")
         ret *= (1 - 2*std::pow(PI*freq*(p[3]-delay),2))*exp(-std::pow(PI*freq*(p[3]-delay),2));
         
      if (decay != 0)
         ret *= std::pow(std::pow(p[0]-center_coords[0],2)+
                     std::pow(p[1]-center_coords[1],2)+
                     std::pow(p[2]-center_coords[2],2),-0.5*decay);
      
      /*Spatial profile of source*/
      if (dir == "rho" || dir == "theta" || dir == "phi")
      {
         Blaze3Vec rtf1({0,0,0}), rtf2({0,0,0});
         if (p[0] != 0 || p[1] != 0 || p[2] != 0)
         {
            rtf1[0] = std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
            rtf1[1] = std::acos(p[2]/rtf1[0]);
            rtf1[2] = std::atan2(p[1],p[0]);
         }
         
         if (center_coords[0] != 0 || center_coords[1] != 0 || center_coords[2] != 0)
         {
            rtf2[0] = std::sqrt(center_coords[0]*center_coords[0]+center_coords[1]*center_coords[1]+center_coords[2]*center_coords[2]);
            rtf2[1] = std::acos(center_coords[2]/rtf2[0]);
            rtf2[2] = std::atan2(center_coords[1],center_coords[0]);
         }
         
         if (prof == "gaussian")
         {
            double exponent = - ((std::pow(p[0]-center_coords[0],2)+
                             std::pow(p[1]-center_coords[1],2)+
                             std::pow(p[2]-center_coords[2],2))/std::pow(radius,2));
            ret *= exp(exponent);
         }
         else if (prof == "wave")
         {
            gsl_sf_result res;
            //~ double nu;
            int lambda;//, status;
            for (uint8_t j=0; j<3; j++)
            {
               switch(bfuncs[j][0])
               {
                  case 'J' : //cylindrical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_Jn_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'Y' : //modified cylindrical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_Yn_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'j' : //spherical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_jl_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'y' : //modified spherical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_yl_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'P' : //Legendre polynomials
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_legendre_Pl_e(lambda,cos(kvec[j]*(rtf1[j]-rtf2[j])),&res);
                     ret *=res.val;
                     break;
                  }
                  default  : 
                  {
                     if (bfuncs[j] == "sin")
                        ret *= sin(2*PI*kvec[j]*(rtf1[j]-rtf2[j]));
                     else if (bfuncs[j] == "cos")
                     {
                        ret *= cos(2*PI*kvec[j]*(rtf1[j]-rtf2[j]));
                     }
                     
                     break;
                  }
               }
            }
         }
         
         if (dir == "rho")
            return Blaze3Vec({sin(rtf1[1])*cos(rtf1[2])*ret,sin(rtf1[1])*sin(rtf1[2])*ret,cos(rtf1[1])*ret});
         else if (dir == "theta")
            return Blaze3Vec({ rtf1[0]*cos(rtf1[1])*cos(rtf1[2])*ret,rtf1[0]*cos(rtf1[1])*sin(rtf1[2])*ret,-rtf1[0]*sin(rtf1[2])*ret});
         else if (dir == "phi")
            return Blaze3Vec({-rtf1[0]*sin(rtf1[1])*sin(rtf1[2])*ret,rtf1[0]*sin(rtf1[1])*cos(rtf1[2])*ret,0});
      }
      else if (dir == "r" || dir == "psi")
      {
         Blaze3Vec rtf1({0,0,0}), rtf2({0,0,0});
         if (p[0] != 0 || p[1] != 0)
         {
            rtf1[0] = std::sqrt(p[0]*p[0]+p[1]*p[1]);
            if (p[0]==0 && p[1] == 0)
               rtf1[1] = 0;
            else if (p[0]>=0)
               rtf1[1] =  std::asin(p[1]/rtf1[0]);
            else
               rtf1[1] = -std::asin(p[1]/rtf1[0])+PI;
            
         }
         
         rtf1[2] = p[2];
         
         if (center_coords[0] != 0 || center_coords[1] != 0 || center_coords[2] != 0)
         {
            rtf2[0] = std::sqrt(center_coords[0]*center_coords[0]+center_coords[1]*center_coords[1]);
            if (center_coords[0]==0 && center_coords[1] == 0)
               rtf2[1] = 0;
            else if (p[0]>=0)
               rtf2[1] =  std::asin(center_coords[1]/rtf2[0]);
            else
               rtf2[1] = -std::asin(center_coords[1]/rtf2[0])+PI;
            // std::cout << "theta = " << rtf1[1]/PI*180 << std::endl;
            
         }
         
         rtf2[2] = center_coords[2];
         
         if (prof == "gaussian")
         {
            double exponent = - ((std::pow(p[0]-center_coords[0],2)+
                             std::pow(p[1]-center_coords[1],2)+
                             std::pow(p[2]-center_coords[2],2))/std::pow(radius,2));
            ret *= exp(exponent);
         }
         else if (prof == "wave")
         {
            // std::cout << "Fai la cosa giusta" << std::endl;
            gsl_sf_result res;
            //~ double nu;
            int lambda;//, status;
            for (uint8_t j=0; j<3; j++)
            {
               switch(bfuncs[j][0])
               {
                  case 'J' : //cylindrical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_Jn_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'Y' : //modified cylindrical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_Yn_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'j' : //spherical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_jl_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'y' : //modified spherical bessel f.
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_bessel_yl_e(lambda,kvec[j]*(rtf1[j]-rtf2[j]),&res);
                     ret *=res.val;
                     break;
                  }
                  case 'P' : //Legendre polynomials
                  {
                     lambda = std::stoi(tail(bfuncs[j],bfuncs[j].size()-1));
                     /*status =*/gsl_sf_legendre_Pl_e(lambda,cos(kvec[j]*(rtf1[j]-rtf2[j])),&res);
                     ret *=res.val;
                     break;
                  }
                  default  : 
                  {
                     if (bfuncs[j] == "sin")
                        ret *= sin(2*PI*kvec[j]*(rtf1[j]-rtf2[j]));
                     else if (bfuncs[j] == "cos")
                     {
                        ret *= cos(2*PI*kvec[j]*(rtf1[j]-rtf2[j]));
                     }
                     
                     break;
                  }
               }
            }
         }
         
         if (dir == "r")
         {
            // std::cout << "----------------------------------------------------------" << std::endl;
            // std::cout << p[0] << '\t' << p[1] << '\t' << p[2] << std::endl;
            // std::cout << rtf1[0] << '\t' << 180*rtf1[1]/PI << '\t' << rtf1[2] << std::endl;
            return Blaze3Vec({cos(rtf1[1])*ret, sin(rtf1[1])*ret,0});
         }
         else if (dir == "psi")
            return Blaze3Vec({ -rtf1[0]*sin(rtf1[1])*ret, rtf1[0]*cos(rtf1[1])*ret, 0});
         else if (dir == "z")
            return Blaze3Vec({0,0,ret});
      }
      else
      {
         if (prof == "gaussian")
         {
            double exponent = - ((std::pow(p[0]-center_coords[0],2)+
                             std::pow(p[1]-center_coords[1],2)+
                             std::pow(p[2]-center_coords[2],2))/std::pow(radius,2));
            ret *= exp(exponent);
         }
         else if (prof == "wave")
         {
            // std::cout << "Fai la cosa giusta" << std::endl;
            
            for (uint8_t j=0; j<3; j++)
            {
               if (bfuncs[j] == "sin")
                  ret *= sin(2*PI*kvec[j]*(p[j]-center_coords[j]));
               else if (bfuncs[j] == "cos")
                  ret *= cos(2*PI*kvec[j]*(p[j]-center_coords[j]));
            }
         }
         
         /*Can handle one principal direction of R^3 per source!
           For arbitrarily oriented sources, 
           use the principle of superposition, i.e. define another source...*/
         if (dir == "x")
            return Blaze3Vec({ret,0,0});
         else if (dir == "y")
            return Blaze3Vec({0,ret,0});
         else
            return Blaze3Vec({0,0,ret});
      }
   }
   
   private:
   Sourcetype st;
   Profile prof;
   Carrier carrier;
   Frequency freq;
   Amplitude amp;
   double    radius, decay, delay, width, impedance;
   Direction dir;
   WaveVector kvec, center_coords;
   Blaze3Vec location;
   std::array<BaseFunction,3> bfuncs;
   std::vector<uint32_t> surface_label;
   bool is_set;
   // double (*foo)();
};
