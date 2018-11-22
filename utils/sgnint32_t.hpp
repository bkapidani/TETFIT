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

// file sgnint32_t<T>.hpp
#ifndef SGNINT32_T_HPP
#define SGNINT32_T_HPP

#include <cmath>
#include <vector>
#include <iostream>
#include <stdlib.h>

template <typename T>
class sgnint32_t
{
   public:
      uint32_t val;
      T sgn;
      sgnint32_t<T>(void) { this->val=0; this->sgn=1; }
      sgnint32_t<T>(uint32_t v, T s) { this->val=v; this->sgn=s; }
      uint32_t Val() const { return this->val; }
      T Sgn()  const { return this->sgn; } 
      sgnint32_t<T>& operator=(const sgnint32_t<T>& i);
      sgnint32_t<T>& operator-(const sgnint32_t<T>& i);
};

// sgnint32_t<T> operator-(const sgnint32_t<T>& i);
// bool operator==(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator!=(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator>(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator<(const sgnint32_t<T> a, const sgnint32_t<T> b);
// bool operator<(const int b, const sgnint32_t<T> a);
// bool operator<(const sgnint32_t<T> a, const int b);
// bool operator==(const sgnint32_t<T> a, const int b);
// bool operator!=(const sgnint32_t<T> a, const int b);
// bool operator<(const sgnint32_t<T> a, const double b);
// uint32_t abs(const sgnint32_t<T> a);
// bool signbit(const sgnint32_t<T>& a);
// int32_t sgn(const sgnint32_t<T>& a);

template <typename T>
sgnint32_t<T>& sgnint32_t<T> :: operator=(const sgnint32_t<T>& i)
{
   this->val=i.Val();
   this->sgn=i.Sgn();
   
   return *this;
}

template <typename T>
sgnint32_t<T> operator-(const sgnint32_t<T>& i)
{
   sgnint32_t<T> ret(i.Val(),-i.Sgn());
   
   return ret;
}

template <typename T>
bool operator==(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval == bval )
      return true;
   else
      return false;
}

template <typename T>
bool operator!=(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval == bval )
      return false;
   else
      return true;
}

template <typename T>
bool operator>(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval > bval )
      return true;
   else
      return false;
}

template <typename T>
bool operator<(const sgnint32_t<T> a, const sgnint32_t<T> b)
{
   auto aval= a.Val();
   auto bval= b.Val();
   
   if ( aval < bval )
      return true;
   else
      return false;
}

template <typename T>
bool operator==(const sgnint32_t<T> a, const int b)
{
   auto aval= a.Val();
   
   if (aval == abs(b))
      return true;
   else
      return false;
}

template <typename T>
bool operator!=(const sgnint32_t<T> a, const int b)
{
   auto aval= a.Val();
   
   if (aval == abs(b))
      return false;
   else
      return true;
}

template <typename T>
bool operator<(const sgnint32_t<T> a, const int b)
{
   auto aval= a.Val();
   auto asgn= a.Sgn();
   
   if ( b == 0 && asgn < 0)
      return true;
   else if (asgn < 0)
      return ((-aval)<b);
   else
      return (aval<b);
}

template <typename T>
bool operator<(const int b, const sgnint32_t<T> a )
{
   auto aval= a.Val();
   auto asgn= a.Sgn();
   
   if ( b == 0 && asgn < 0)
      return false;
   else if (asgn < 0)
      return (b<(-aval));
   else
      return (b<aval);
}

template <typename T>
bool operator<(const sgnint32_t<T> a, const double b)
{
   auto aval= a.Val();
   auto asgn= a.Sgn();
   
   if ( b == 0 && asgn < 0)
      return true;
   else if (asgn < 0)
      return ((-aval)<b);
   else
      return (aval<b);
}

template <typename T>
uint32_t abs(const sgnint32_t<T> a)
{   
   return a.Val(); 
}

template <typename T>
bool signbit(const sgnint32_t<T>& a)
{
   return a.Sgn()<0; 
}

template <typename T>
T sgn(const sgnint32_t<T>& a)
{
   return a.Sgn();
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const sgnint32_t<T>& s)
{
    T d = s.Sgn()*(s.Val());
    os << d;
    
    return os;
}

template <typename T>
std::vector<uint32_t> abs(const std::vector<sgnint32_t<T>> v)
{
   std::vector<uint32_t> ret;
   for (auto vv : v)
      ret.push_back(abs(vv));
   return ret;
      
} 

#endif
