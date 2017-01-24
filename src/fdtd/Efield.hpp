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
 
#ifndef EFIELD_HPP
#define EFIELD_HPP

#include "timecounter.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

template<typename T>
class Hfield;

template<typename T>
class Efield
{
   // friend class H_field;
   public:
   Efield(double t, int32_t c, const std::vector<uint32_t> C, bool is_b, T bc)
    : t_step(t), curl(c), C(C), is_boundary(is_b), bc(bc)
   {
	   value = 0;
   }
   
   T& Update(const std::vector<Hfield<T>>& F, const T& time_function)
   {
	   value = is_boundary ? time_function*bc : value + t_step*M*curl*(F[C[0]]-F[C[1]]+F[C[2]]-F[C[3]]);
	   return value;
   }
   
   const std::vector<uint32_t> C;
   T value, M, bc;
   int32_t curl;
   double t_step;
   bool is_boundary;
};

#endif
