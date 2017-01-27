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
 
#ifndef MESHER_HPP
#define MESHER_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "timecounter.h"
#include "Efield.hpp"
#include "Hfield.hpp"
#include <silo.h>
#include <string.h>
#include <set>
#include <cassert>
#include <fstream>
#include <iomanip>

/*general parameters*/
double pi = 3.141592653589793;
double mu0 = 4*pi*1e-7;
double epsilon0 = 8.854187817e-12;
double c0 = 1 / sqrt( mu0 * epsilon0 );
Eigen::Vector3d probe_point(0.025,0.0125,0.05);

template<typename T>
class mesher
{
   public:
   mesher(void)
   {
	   timecounter t_mesh;
	   t_mesh.tic();
	   uint32_t nv,nf,ne,np;
	   nv=nf=ne=np=0;
	   epsilon[1]=epsilon0;
	   mu[1]=mu0;
	   freq=5e9;
	   
      //Numerical limits
		 
      xmin= 0;
      xmax= 0.05;
      ymin= 0;
      ymax= 0.025;
      zmin= 0;
      zmax= 0.1;
      L=0.0005;
      //other instructions in constructor
	  // std::cout << "Qui?" << std::endl;
	  auto px = xmin;
	  auto py = ymin;
	  auto pz = zmin;  
	  
	  Lx = Ly = Lz = L;
	  t_step = 0.5*sqrt(pow(Lx,2)+pow(Ly,2)+pow(Lz,2))/c0/sqrt(3);
	  volume = Lx*Ly*Lz;

	  T area_z = Lx*Ly;
	  T area_x = Ly*Lz;
	  T area_y = Lx*Lz;
	  
	  
	  
	  T da_z = area_z/4;
	  T da_x = area_x/4;
	  T da_y = area_y/4;

	  Eigen::Vector3d area_z_vec(0,0,Lx*Ly);
	  Eigen::Vector3d area_y_vec(0,Lx*Lz,0);
	  Eigen::Vector3d area_x_vec(Ly*Lz,0,0);
	  
	  Eigen::Vector3d dual_area_z(0,0,Lx*Ly/4);
	  Eigen::Vector3d dual_area_y(0,Lx*Lz/4,0);
	  Eigen::Vector3d dual_area_x(Ly*Lz/4,0,0);
	  
	  this->dual_area_x=dual_area_x;
	  this->dual_area_y=dual_area_y;
	  this->dual_area_z=dual_area_z;

	  this->area_x_vec=area_x_vec;
	  this->area_y_vec=area_y_vec;
	  this->area_z_vec=area_z_vec;
	  
	  Nx = (fabs(xmax-xmin)) / Lx;// + 1;
	  Ny = (fabs(ymax-ymin)) / Ly;// + 1;
	  Nz = (fabs(zmax-zmin)) / Lz;// + 1;
	  
	  tot_E= Nx*(Ny+1)*(Nz+1)+(Nx+1)*Ny*(Nz+1)+(Nx+1)*(Ny+1)*Nz;
	  tot_F= Nx*Ny*(Nz+1)+Nz*Nx*(Ny+1)+Ny*Nz*(Nx+1);
	  
	  typedef Eigen::Triplet<uint32_t> Q;
      
	  Eigen::SparseMatrix<uint32_t> previous_layer(Nx,Ny);
	  std::vector<T> average_eps(tot_E,0), average_ni(tot_F,0), face_area(tot_F), edge_len(tot_E);
	  // previous_layer.SetZero();
	
      Eigen::Vector3d inc_x(Lx,0,0), inc_y(0,Ly,0), inc_z(0,0,Lz), dummy_vec;
	  bool not_found = true;
	  std::vector<int32_t> cplus=std::vector<int32_t>({1,-1,1,-1});
	  std::vector<int32_t> cminus=std::vector<int32_t>({-1,1,-1,1});
	  
	  for(uint32_t k=0;k<Nz;k++)
      {
         py=ymin;
		 std::vector<uint32_t> old_col;
		 Eigen::SparseMatrix<uint32_t> this_layer(Nx,Ny);
         std::vector<Q> tripletList;
         tripletList.reserve(Nx*Ny);

         for (uint32_t j=0;j<Ny;j++)
		 {
            std::vector<uint32_t> this_col(Nx,0);
            px=xmin;
            for (uint32_t i=0;i<Nx;i++)
			{
			   
               if (IsGridVol(px,py,pz))
			   {
				  Eigen::Vector3d pp(px,py,pz);
				  
				  if (not_found)
				  {
					 auto pp2 = pp+inc_x+inc_y+inc_z;
				     if (probe_point(0)>= pp(0) && probe_point(0)<= pp2(0) && 
					     probe_point(1)>= pp(1) && probe_point(1)<= pp2(1) &&
                         probe_point(2)>= pp(2) && probe_point(2)<= pp2(2))
					 {
                        not_found = false;
						probe_elem = nv;
						// std::cout << "# of element to be probed: " << nv+1 << std::endl;
					 }
                  }

				  std::vector<uint32_t> dummy(6), dummyf;
				  std::vector<int32_t> dummycurl;
				  uint8_t boxtype=0;
				  uint32_t bottom,left,back;
				  bottom = left = back = 0;
				  
                  if (pz>zmin)
					 bottom = previous_layer.coeffRef(i,j);
				  if (py>ymin)
					 left = old_col[i];
                  if (px>xmin)
					 back = this_col[i-1];
				 
				  if (bottom)
					  boxtype++;
				  if (left)
					  boxtype+=2;
				  if (back)
					  boxtype+=4;
				  
				  switch(boxtype)
				  {
                     case 0 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
                        D.push_back(std::vector<uint32_t>({nf,nf+1,nf+2,
						                                      nf+3,nf+4,nf+5}));
						nf+=6;
						
						for (uint32_t cnt=0; cnt<6; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						
						E_cluster.push_back(std::vector<uint32_t>({ne,ne+1,ne+2,ne+3,ne+4,ne+5,
						                                    ne+6,ne+7,ne+8,ne+9,ne+10,ne+11}));
						ne+=12;
						
						for (uint32_t cnt=0; cnt<12; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						
						P_cluster.push_back(std::vector<uint32_t>({np,np+1,np+2,np+3,
						                                    np+4,np+5,np+6,np+7}));
						np+=8;
						
						// for (uint32_t cnt=0; cnt<8; cnt++)
                           //Gt.push_back(dummyf);
						
                        pts.push_back(pp);
						pts.push_back(pp+inc_x);
						pts.push_back(pp+inc_y);
						pts.push_back(pp+inc_x+inc_y);
                        pts.push_back(pp+inc_z);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 1 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<uint32_t>({D[bottom-1][5],nf,nf+1,
						                                      nf+2,nf+3,nf+4}));
						nf+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],ne,E_cluster[bottom-1][10],
						                                    ne+1,E_cluster[bottom-1][11],ne+2,ne+3,ne+4,ne+5,ne+6,ne+7}));
						ne+=8;
						
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}

						P_cluster.push_back(std::vector<uint32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    np,np+1,np+2,np+3}));
						np+=4;
						
 						// for (uint32_t cnt=0; cnt<4; cnt++)
							//Gt.push_back(dummyf);
						
                        pts.push_back(pp+inc_z);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 2 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<uint32_t>({nf,D[left-1][4],nf+1,
						                                      nf+2,nf+3,nf+4}));
						nf+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({E_cluster[left-1][5],ne,E_cluster[left-1][6],ne+1,E_cluster[left-1][7],ne+2,
						                                    ne+3,ne+4,E_cluster[left-1][11],ne+5,ne+6,ne+7}));
						ne+=8;
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<uint32_t>({P_cluster[left-1][2],P_cluster[left-1][3],np,np+1,
						                                    P_cluster[left-1][6],P_cluster[left-1][7],np+2,np+3}));
						np+=4;
						
						// for (uint32_t cnt=0; cnt<4; cnt++)
							//Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_y);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 3 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
                        D.push_back(std::vector<uint32_t>({D[bottom-1][5],D[left-1][4],nf,
						                                      nf+1,nf+2,nf+3}));
						nf+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],E_cluster[left-1][6],E_cluster[bottom-1][10],
						                                          E_cluster[left-1][7],E_cluster[bottom-1][11],ne,ne+1,
																  E_cluster[left-1][11],ne+2,ne+3,ne+4}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<uint32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                          P_cluster[left-1][6],P_cluster[left-1][7],np,np+1}));
						np+=2;
						
						// for (uint32_t cnt=0; cnt<2; cnt++)
							//Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 4 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<uint32_t>({nf,nf+1,D[back-1][3],
						                                      nf+2,nf+3,nf+4}));
						nf+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({ne,E_cluster[back-1][3],E_cluster[back-1][4],ne+1,ne+2,ne+3,
						                                    E_cluster[back-1][7],ne+4,ne+5,E_cluster[back-1][10],ne+6,ne+7}));
						ne+=8;
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<uint32_t>({P_cluster[back-1][1],np,P_cluster[back-1][3],np+1,
						                                    P_cluster[back-1][5],np+2,P_cluster[back-1][7],np+3}));
						np+=4;
						
						// for (uint32_t cnt=0; cnt<4; cnt++)
							//Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_x);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 5 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<uint32_t>({D[bottom-1][5],nf,D[back-1][3],
						                                      nf+1,nf+2,nf+3}));
						nf+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],E_cluster[back-1][4],E_cluster[bottom-1][10],
						                                          ne,E_cluster[bottom-1][11],E_cluster[back-1][7],ne+1,
																  ne+2,E_cluster[back-1][10],ne+3,ne+4}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<uint32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    P_cluster[back-1][5],np,P_cluster[back-1][7],np+1}));
						np+=2;
						
						// for (uint32_t cnt=0; cnt<2; cnt++)
							//Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 6 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<uint32_t>({nf,D[left-1][4],D[back-1][3],
						                                      nf+1,nf+2,nf+3}));
						nf+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({E_cluster[left-1][5],E_cluster[back-1][3],E_cluster[left-1][6],ne,
						                                          E_cluster[left-1][7],ne+1,E_cluster[back-1][7],ne+2,
																  E_cluster[left-1][11],E_cluster[back-1][10],ne+3,ne+4}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<uint32_t>({P_cluster[left-1][2],P_cluster[left-1][3],P_cluster[back-1][3],np,
						                                    P_cluster[left-1][6],P_cluster[left-1][7],P_cluster[back-1][7],np+1}));
						np+=2;
						
						// for (uint32_t cnt=0; cnt<2; cnt++)
							//Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 7 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
                        D.push_back(std::vector<uint32_t>({D[bottom-1][5],D[left-1][4],D[back-1][3],
						                                      nf,nf+1,nf+2}));
						nf+=3;
						for (uint32_t cnt=0; cnt<3; cnt++)
						{
							//Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							
						}
						E_cluster.push_back(std::vector<uint32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],E_cluster[left-1][6],E_cluster[bottom-1][10],
						                                          E_cluster[left-1][7],E_cluster[bottom-1][11],E_cluster[back-1][7],ne,
																  E_cluster[left-1][11],E_cluster[back-1][10],ne+1,ne+2}));
						ne+=3;
						for (uint32_t cnt=0; cnt<3; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<uint32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    P_cluster[left-1][6],P_cluster[left-1][7],P_cluster[back-1][7],np}));
						np+=1;
						
						//Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
				  }
				  
				  
				  
                  //Dt[D[nv][0]].push_back( nv); 
                  //Dt[D[nv][1]].push_back( nv);
				  //Dt[D[nv][2]].push_back( nv);
				  //Dt[D[nv][3]].push_back( nv);
				  //Dt[D[nv][4]].push_back( nv);
				  //Dt[D[nv][5]].push_back( nv);

				  material.push_back(1);
				
				 if (!G[E_cluster[nv][0]].size())
				 {
					G[E_cluster[nv][0]]  = std::vector<uint32_t>({P_cluster[nv][0],P_cluster[nv][1]});
					edge_len[E_cluster[nv][0]] = Lx;
					//Gt[P_cluster[nv][0]].push_back(E_cluster[nv][0]);
					//Gt[P_cluster[nv][1]].push_back( E_cluster[nv][0]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][1]].size())
				 {
					G[E_cluster[nv][1]]  = std::vector<uint32_t>({P_cluster[nv][0],P_cluster[nv][2]});
					edge_len[E_cluster[nv][1]] = Ly;
					//Gt[P_cluster[nv][0]].push_back(E_cluster[nv][1]);
					//Gt[P_cluster[nv][2]].push_back( E_cluster[nv][1]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][2]].size())
				 {
					G[E_cluster[nv][2]]  = std::vector<uint32_t>({P_cluster[nv][0],P_cluster[nv][4]});
					edge_len[E_cluster[nv][2]] = Lz;
					//Gt[P_cluster[nv][0]].push_back(E_cluster[nv][2]);
					//Gt[P_cluster[nv][4]].push_back( E_cluster[nv][2]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][3]].size())
				 {
					G[E_cluster[nv][3]]  = std::vector<uint32_t>({P_cluster[nv][1],P_cluster[nv][3]});
					edge_len[E_cluster[nv][3]] = Ly;
					//Gt[P_cluster[nv][1]].push_back(E_cluster[nv][3]);
					//Gt[P_cluster[nv][3]].push_back( E_cluster[nv][3]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][4]].size())
				 {
					G[E_cluster[nv][4]]  = std::vector<uint32_t>({P_cluster[nv][1],P_cluster[nv][5]});
					edge_len[E_cluster[nv][4]] = Lz;
					//Gt[P_cluster[nv][1]].push_back(E_cluster[nv][4]);
					//Gt[P_cluster[nv][5]].push_back( E_cluster[nv][4]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][5]].size())
				 {
					G[E_cluster[nv][5]]  = std::vector<uint32_t>({P_cluster[nv][2],P_cluster[nv][3]});
					edge_len[E_cluster[nv][5]] = Lx;
					//Gt[P_cluster[nv][2]].push_back(E_cluster[nv][5]);
					//Gt[P_cluster[nv][3]].push_back( E_cluster[nv][5]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][6]].size())
				 {
					G[E_cluster[nv][6]]  = std::vector<uint32_t>({P_cluster[nv][2],P_cluster[nv][6]});
					edge_len[E_cluster[nv][6]] = Lz;
					//Gt[P_cluster[nv][2]].push_back(E_cluster[nv][6]);
					//Gt[P_cluster[nv][6]].push_back( E_cluster[nv][6]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][7]].size())
				 {
					G[E_cluster[nv][7]]  = std::vector<uint32_t>({P_cluster[nv][3],P_cluster[nv][7]});
					edge_len[E_cluster[nv][7]] = Lz;
					//Gt[P_cluster[nv][3]].push_back(E_cluster[nv][7]);
					//Gt[P_cluster[nv][7]].push_back( E_cluster[nv][7]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][8]].size())
				 {
					G[E_cluster[nv][8]]  = std::vector<uint32_t>({P_cluster[nv][4],P_cluster[nv][5]});
					edge_len[E_cluster[nv][8]] = Lx;
					//Gt[P_cluster[nv][4]].push_back(E_cluster[nv][8]);
					//Gt[P_cluster[nv][5]].push_back( E_cluster[nv][8]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][9]].size())
				 {
					G[E_cluster[nv][9]]  = std::vector<uint32_t>({P_cluster[nv][4],P_cluster[nv][6]});
					edge_len[E_cluster[nv][9]] = Ly;
					//Gt[P_cluster[nv][4]].push_back(E_cluster[nv][9]);
					//Gt[P_cluster[nv][6]].push_back( E_cluster[nv][9]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][10]].size())
				 {
					G[E_cluster[nv][10]]  = std::vector<uint32_t>({P_cluster[nv][5],P_cluster[nv][7]});
					edge_len[E_cluster[nv][10]] = Ly;
					//Gt[P_cluster[nv][5]].push_back(E_cluster[nv][10]);
					//Gt[P_cluster[nv][7]].push_back( E_cluster[nv][10]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv][11]].size())
				 {
					G[E_cluster[nv][11]]  = std::vector<uint32_t>({P_cluster[nv][6],P_cluster[nv][7]});
					edge_len[E_cluster[nv][11]] = Lx;
					//Gt[P_cluster[nv][6]].push_back(E_cluster[nv][11]);
					//Gt[P_cluster[nv][7]].push_back( E_cluster[nv][11]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 
				 if (!C[D[nv][0]].size())
				 {
					C[D[nv][0]] = std::vector<uint32_t>({E_cluster[nv][0],E_cluster[nv][1],E_cluster[nv][3],E_cluster[nv][5]});
					
					Ct[E_cluster[nv][0]].push_back(D[nv][0]);
					Ct[E_cluster[nv][1]].push_back(D[nv][0]);
					Ct[E_cluster[nv][3]].push_back(D[nv][0]);
					Ct[E_cluster[nv][5]].push_back(D[nv][0]);

					face_area[D[nv][0]] = area_z;
					F.push_back(0);
					boundary_face.push_back(1);
				 }
				 else
					 boundary_face[D[nv][0]]=0;
					
				 
				 if (!C[D[nv][1]].size())
				 {
					C[D[nv][1]] = std::vector<uint32_t>({E_cluster[nv][0],E_cluster[nv][2],E_cluster[nv][4],E_cluster[nv][8]});

					Ct[E_cluster[nv][0]].push_back(D[nv][1]);
					Ct[E_cluster[nv][2]].push_back(D[nv][1]);
					Ct[E_cluster[nv][4]].push_back(D[nv][1]);
					Ct[E_cluster[nv][8]].push_back(D[nv][1]);

					curl[D[nv][1]] = -1;
					face_area[D[nv][1]] = area_y;
					F.push_back(0);
					boundary_face.push_back(1);
				 }
				 else
					 boundary_face[D[nv][1]]=0;
				 
				 if (!C[D[nv][2]].size())
				 {
					C[D[nv][2]] = std::vector<uint32_t>({E_cluster[nv][1],E_cluster[nv][2],E_cluster[nv][6],E_cluster[nv][9]});

					Ct[E_cluster[nv][1]].push_back(D[nv][2]);
					Ct[E_cluster[nv][2]].push_back(D[nv][2]);
					Ct[E_cluster[nv][6]].push_back(D[nv][2]);
					Ct[E_cluster[nv][9]].push_back(D[nv][2]);

					face_area[D[nv][2]] = area_x;					
					F.push_back(0);
					boundary_face.push_back(1);
				 }
				 else
					 boundary_face[D[nv][2]]=0;
				 
				 if (!C[D[nv][3]].size())
				 {
					C[D[nv][3]] = std::vector<uint32_t>({E_cluster[nv][3],E_cluster[nv][4],E_cluster[nv][7],E_cluster[nv][10]});

					Ct[E_cluster[nv][3]].push_back(D[nv][3]);
					Ct[E_cluster[nv][4]].push_back(D[nv][3]);
					Ct[E_cluster[nv][7]].push_back(D[nv][3]);
					Ct[E_cluster[nv][10]].push_back(D[nv][3]);

					face_area[D[nv][3]] = area_x;					
					F.push_back(0);
					boundary_face.push_back(1);
				 }
				 else
					 boundary_face[D[nv][3]]=0;
				 
				 if (!C[D[nv][4]].size())
				 {
					C[D[nv][4]] = std::vector<uint32_t>({E_cluster[nv][5],E_cluster[nv][6],E_cluster[nv][7],E_cluster[nv][11]});

					Ct[E_cluster[nv][ 5]].push_back(D[nv][4]);
					Ct[E_cluster[nv][ 6]].push_back(D[nv][4]);
					Ct[E_cluster[nv][ 7]].push_back(D[nv][4]);
					Ct[E_cluster[nv][11]].push_back(D[nv][4]);

					curl[D[nv][4]] = -1;
					face_area[D[nv][4]] = area_y;					
					F.push_back(0);
					boundary_face.push_back(1);
				 }
				 else
					 boundary_face[D[nv][4]]=0;
				 
				 if (!C[D[nv][5]].size())
				 {
					C[D[nv][5]] = std::vector<uint32_t>({E_cluster[nv][8],E_cluster[nv][9],E_cluster[nv][10],E_cluster[nv][11]});

					Ct[E_cluster[nv][ 8]].push_back(D[nv][5]); 
					Ct[E_cluster[nv][ 9]].push_back(D[nv][5]);
					Ct[E_cluster[nv][10]].push_back(D[nv][5]);
					Ct[E_cluster[nv][11]].push_back(D[nv][5]);

					face_area[D[nv][5]] = area_z;					
					F.push_back(0);
					boundary_face.push_back(1);
				 }
				 else
					 boundary_face[D[nv][5]]=0;

					average_ni[D[nv][0]] += Lz/2/mu[material[nv]];///face_area[D[nv][0]];
					average_ni[D[nv][1]] += Ly/2/mu[material[nv]];///face_area[D[nv][1]];
					average_ni[D[nv][2]] += Lx/2/mu[material[nv]];///face_area[D[nv][2]];
					average_ni[D[nv][3]] += Lx/2/mu[material[nv]];///face_area[D[nv][3]];
					average_ni[D[nv][4]] += Ly/2/mu[material[nv]];///face_area[D[nv][4]];
					average_ni[D[nv][5]] += Lz/2/mu[material[nv]];///face_area[D[nv][5]];
					
					average_eps[E_cluster[nv][ 0]] += da_x*epsilon[material[nv]];///edge_len[E_cluster[nv][ 0]];
					average_eps[E_cluster[nv][ 1]] += da_y*epsilon[material[nv]];///edge_len[E_cluster[nv][ 1]];
					average_eps[E_cluster[nv][ 2]] += da_z*epsilon[material[nv]];///edge_len[E_cluster[nv][ 2]];
					average_eps[E_cluster[nv][ 3]] += da_y*epsilon[material[nv]];///edge_len[E_cluster[nv][ 3]];
					average_eps[E_cluster[nv][ 4]] += da_z*epsilon[material[nv]];///edge_len[E_cluster[nv][ 4]];
					average_eps[E_cluster[nv][ 5]] += da_x*epsilon[material[nv]];///edge_len[E_cluster[nv][ 5]];
					average_eps[E_cluster[nv][ 6]] += da_z*epsilon[material[nv]];///edge_len[E_cluster[nv][ 6]];
					average_eps[E_cluster[nv][ 7]] += da_z*epsilon[material[nv]];///edge_len[E_cluster[nv][ 7]];
					average_eps[E_cluster[nv][ 8]] += da_x*epsilon[material[nv]];///edge_len[E_cluster[nv][ 8]];
					average_eps[E_cluster[nv][ 9]] += da_y*epsilon[material[nv]];///edge_len[E_cluster[nv][ 9]];
					average_eps[E_cluster[nv][10]] += da_y*epsilon[material[nv]];///edge_len[E_cluster[nv][10]];
					average_eps[E_cluster[nv][11]] += da_x*epsilon[material[nv]];///edge_len[E_cluster[nv][11]];
					
					// dual_pts.push_back(pp+0.5*inc_x+0.5*inc_y+0.5*inc_z);
					nv++;
					tripletList.push_back(Q(i,j,nv));
					this_col[i]=nv;
			   }
			   
               px+=Lx;
			}
			
			old_col=std::move(this_col);
			py+=Ly;
		 }

         this_layer.setFromTriplets(tripletList.begin(), tripletList.end());
		 previous_layer=std::move(this_layer);
		 pz+=Lz;
      }
	  
	  for (uint32_t i=0; i<nf; i++)
		M_ni.push_back(average_ni[i]/face_area[i]);
	  
	  for (uint32_t i=0; i<ne; i++)
	  {
		 M_h.push_back(edge_len[i]/average_eps[i]);
		 uint8_t in_b=0;
		 for (auto ff : Ct[i])
		 {
			if (boundary_face[ff])
			{
				in_b = 1;
				this->is_boundary.push_back(true);
				bc_edges.push_back(i);
				break;
			}
		 }
		 
		 switch (in_b)
		 {
			case  1 :
			{
				// this->is_boundary.push_back(true);

				auto bar = edge_barycenter(i);
				if (bar(2) <= zmin && Ct[i].size() > 2)
					bc.push_back(sin( bar(0)*pi/0.05 )*edge_vector(i)(1));
				else
					bc.push_back(0);

				break;
			}
			case 0 :
			{
				bc.push_back(1);
				this->is_boundary.push_back(false);

				break;
			}
			default  :
			{
				std::cout << "Wrong adjacency list routine!" << std::endl;
				break;
			}
		 }
	  }
	  
      t_mesh.toc();
	  std::cout << "Meshing and material modeling done in " << t_mesh << " seconds" << std::endl;
	  std::cout << "Mesh statistics: " << std::endl;
	  std::cout << std::setw(20) << "volumes: "           << std::setw(10) << nv << std::endl; 
	  std::cout << std::setw(20) << "surfaces: "          << std::setw(10) << nf << std::endl;
	  std::cout << std::setw(20) << "edges: "             << std::setw(10) << ne << std::endl; 
	  std::cout << std::setw(20) << "vertices: "          << std::setw(10) << np << std::endl;
	  std::cout << std::setw(20) << "materials: "         << std::setw(10) << mu.size() << std::endl; 
	  std::cout << std::setw(20) << "N. of unknowns: "    << std::setw(10) << nf+ne << std::endl;
   }
   
   void Run_VoxBased(const double simulation_time)
   {
	  std::cout << "---------------- Running FDTD simulation ----------------" << std::endl << std::endl;
	  double step_time_average=0;
	  const uint32_t N_of_steps=simulation_time/t_step;
	  size_t i;
	  
	  std::cout << Nx << " " << Ny << " " << Nz << std::endl;
	  
	  T time_function;
	  timecounter step_cost;
	  std::vector<T> numeric_values,numeric_times;
	  
		for (i=0; i*t_step <= simulation_time; i++)
		{
			step_cost.tic();
			
			/* Handle boundary field excitations */
			time_function=sin(2*pi*freq*i*t_step);
			for (auto ee : bc_edges)
				if (bc[ee])
					U[ee] = time_function*bc[ee];
			
			
			/* Handle fields voxel per voxel */
			for (uint32_t j = 0; j<D.size()-1; j++)
			{
				if (!is_boundary[E_cluster[j][ 7]])
					U[E_cluster[j][ 7]] += t_step*M_h[E_cluster[j][ 7]]*(F[D[j][3]]-F[D[j+Nx][3]]+F[D[j+1][4]]-F[D[j][4]]);
				if (!is_boundary[E_cluster[j][10]])
					U[E_cluster[j][10]] += t_step*M_h[E_cluster[j][10]]*(F[D[j][5]]-F[D[j+1][5]]+F[D[j+Nx*Ny][3]]-F[D[j][3]]);
				if (!is_boundary[E_cluster[j][11]])
					U[E_cluster[j][11]] += t_step*M_h[E_cluster[j][11]]*(F[D[j][4]]-F[D[j+Nx*Ny][4]]+F[D[j+Nx][5]]-F[D[j][5]]);

				if (!boundary_face[D[j][3]])
					F[D[j][3]] -= t_step*M_ni[D[j][3]]*(U[E_cluster[j][7]]-U[E_cluster[j][4]]+U[E_cluster[j][3]]-U[E_cluster[j][10]]);
				if (!boundary_face[D[j][4]])
					F[D[j][4]] -= t_step*M_ni[D[j][4]]*(U[E_cluster[j][6]]-U[E_cluster[j][7]]+U[E_cluster[j][11]]-U[E_cluster[j][5]]);
				if (!boundary_face[D[j][5]])
					F[D[j][5]] -= t_step*M_ni[D[j][5]]*(U[E_cluster[j][8]]-U[E_cluster[j][11]]+U[E_cluster[j][10]]-U[E_cluster[j][9]]);
			}
			
			auto num_val = GetElectricField(probe_elem);
			numeric_values.push_back(num_val(1));
			numeric_times.push_back(i*t_step);
			step_cost.toc();
			step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
			
			if (i % 20 == 0)
				ExportFields(i);
			
			if ((i+1) % 140 == 0)
				std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(7) << step_time_average << "s, " 
			              << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
		}
		
		/* Output stats and fields*/
		std::ofstream os;
		os.open("numeric_FIT.dat");
		for (size_t k=0; k < numeric_values.size(); k++)
		  os << numeric_times[k] << " " << numeric_values[k] << std::endl;
		os.close();
		std::cout << "Time step takes (average) " << step_time_average/(double(i)) << " seconds (" << i << " time steps!)" << std::endl;
		std::cout << "Total running time is "     << step_time_average << " seconds" << std::endl;
   }
   
   void Run(const double simulation_time)
   {
	  std::cout << "---------------- Running FDTD simulation ----------------" << std::endl << std::endl;
	  double step_time_average=0;
	  const uint32_t N_of_steps=simulation_time/t_step;
	  size_t i;
	  
	  T time_function;
	  timecounter step_cost;
	  std::vector<T> numeric_values,numeric_times;

      for (i=0; i*t_step <= simulation_time; i++)
	  {
		 step_cost.tic();
		 time_function=sin(2*pi*freq*i*t_step);
         for (size_t j=0; j<U.size(); j++)
		 {
			if (bc[j]!=0)
			{
				if (!is_boundary[j])
					U[j] += t_step*M_h[j]*dual_curl[j]*(F[Ct[j][0]]-F[Ct[j][1]]+F[Ct[j][2]]-F[Ct[j][3]]);
				else
					U[j] = time_function*bc[j];
			}
		 }
		 
         for (size_t j=0; j<F.size(); j++)
			 if (!boundary_face[j])
				F[j] -= t_step*M_ni[j]*curl[j]*(U[C[j][0]]-U[C[j][1]]+U[C[j][2]]-U[C[j][3]]);
		
		 auto num_val = GetElectricField(probe_elem);
		 numeric_values.push_back(num_val(1));
		 numeric_times.push_back(i*t_step);
		 step_cost.toc();
		 step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
		 
		if (i % 20 == 0)
			ExportFields(i);
		
		if ((i+1) % 140 == 0)
			std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(7) << step_time_average << "s, " 
					  << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
	  }
	  
	  /* Output stats and fields*/
	  std::ofstream os;
	  os.open("numeric_FIT.dat");
	  for (size_t k=0; k < numeric_values.size(); k++)
		  os << numeric_times[k] << " " << numeric_values[k] << std::endl;
	  os.close();
	  std::cout << "Time step takes (average) " << step_time_average/(double(i)) << " seconds (" << i << " time steps!)" << std::endl;
	  std::cout << "Total running time is "     << step_time_average << " seconds" << std::endl;
   }
   
   uint32_t Volumes_size() { return D.size(); }
   uint32_t Surfaces_size() { return C.size(); }
   uint32_t Edges_size() { return G.size(); }
   uint32_t Points_size() { return pts.size(); }

	void ExportFields (uint32_t time)
	{
		timecounter t_export;
		t_export.tic();
		std::stringstream namedummy;
		namedummy << "./output/fit" << time << ".silo";
		auto filename = namedummy.str();
		_siloDb = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
		
		auto Meshname = ExportMesh(filename);
		
		std::vector<double> Ex_vals, Ey_vals, Ez_vals, Hx_vals, Hy_vals, Hz_vals;
		Ex_vals.reserve( Volumes_size() );
		Ey_vals.reserve( Volumes_size() );
		Ez_vals.reserve( Volumes_size() );
		Hx_vals.reserve( Volumes_size() );
		Hy_vals.reserve( Volumes_size() );
		Hz_vals.reserve( Volumes_size() );
		
		for (auto itor = 0; itor < Volumes_size(); itor++)
		{
			auto Efield = GetElectricField(itor);
			auto Hfield = GetMagneticField(itor);
			Ex_vals.push_back(Efield(0));
			Ey_vals.push_back(Efield(1));
			Ez_vals.push_back(Efield(2));
			Hx_vals.push_back(Efield(0));
			Hy_vals.push_back(Efield(1));
			Hz_vals.push_back(Efield(2));
		}
		
		int dims[3] = { static_cast<int>(Nx), static_cast<int>(Ny), static_cast<int>(Nz) };
		
		std::vector<std::string> varnames({"Ex","Ey","Ez","Hx","Hy","Hz"});
		DBPutQuadvar1(_siloDb, varnames[0].c_str(), Meshname.c_str(), Ex_vals.data(),dims,3,NULL,0,DB_DOUBLE, DB_ZONECENT, NULL);
		DBPutQuadvar1(_siloDb, varnames[1].c_str(), Meshname.c_str(), Ey_vals.data(),dims,3,NULL,0,DB_DOUBLE, DB_ZONECENT, NULL);
		DBPutQuadvar1(_siloDb, varnames[2].c_str(), Meshname.c_str(), Ez_vals.data(),dims,3,NULL,0,DB_DOUBLE, DB_ZONECENT, NULL);
		DBPutQuadvar1(_siloDb, varnames[3].c_str(), Meshname.c_str(), Hx_vals.data(),dims,3,NULL,0,DB_DOUBLE, DB_ZONECENT, NULL);
		DBPutQuadvar1(_siloDb, varnames[4].c_str(), Meshname.c_str(), Hy_vals.data(),dims,3,NULL,0,DB_DOUBLE, DB_ZONECENT, NULL);
		DBPutQuadvar1(_siloDb, varnames[5].c_str(), Meshname.c_str(), Hz_vals.data(),dims,3,NULL,0,DB_DOUBLE, DB_ZONECENT, NULL);
		
		DBClose(_siloDb);
		
		t_export.toc();
		std::cout << "SILO: Output to file done in " << std::setw(7) << t_export << std::setw(8) << " seconds" << std::endl;
	}

	std::string ExportMesh(std::string filename)
	{
		// DBfile *file = NULL; /* The Silo file pointer */
		int32_t tot_pt = Points_size();
		
		std::stringstream meshname;
		
		meshname << "fit_grid_" << tot_pt << "nodes";
		
		char *coordnames[3]; /* Names of the coordinates */
		std::vector<double> nodex, nodey, nodez;
		nodex.reserve(tot_pt);
		nodey.reserve(tot_pt);
		nodez.reserve(tot_pt);
		
		/* Name the coordinate axes ‘X’ and ‘Y’ */
		coordnames[0] = strdup("X");
		coordnames[1] = strdup("Y");
		coordnames[2] = strdup("Z");
		/* Give the x coordinates of the mesh */
		// uint32_t pt_iter=0;
		for (size_t i=0; i<=Nx; i++)
			nodex.push_back(i*Lx);
		for (size_t i=0; i<=Ny; i++)
			nodey.push_back(i*Ly);
		for (size_t i=0; i<=Nz; i++)
			nodez.push_back(i*Lz);
		
		/* How many nodes in each direction? */
		int dimensions[3]={static_cast<int>(Nx+1),static_cast<int>(Ny+1),static_cast<int>(Nz+1)};
		/* Assign coordinates to coordinates array */
		double *coordinates[3];
		coordinates[0] = nodex.data();
		coordinates[1] = nodey.data();
		coordinates[2] = nodez.data();
		
        std::vector<int> nodelist;
        nodelist.reserve( Volumes_size() );
        
		DBPutQuadmesh(_siloDb,meshname.str().c_str(), coordnames,coordinates,dimensions,3,DB_DOUBLE, DB_COLLINEAR, NULL);

		return meshname.str();
	}
   
   private:
   T xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz,L,volume;
   uint32_t probe_elem, tot_E, tot_F, Nx, Ny, Nz;
   std::vector<bool> is_boundary;
   std::vector<T> M_ni, M_h, U, F, bc;
   std::map<uint32_t,double> epsilon,mu;
   double t_step, freq;
   DBfile *_siloDb=NULL;
   std::vector<uint32_t> material, bc_edges;
   std::vector<uint8_t> boundary_face;
   std::vector<std::vector<uint32_t>> D/*,Dt*/,C,Ct,G/*,Gt*/;
   std::vector<int32_t> curl, dual_curl;
   Eigen::Vector3d dual_area_z, dual_area_y, dual_area_x;
   Eigen::Vector3d area_z_vec, area_y_vec, area_x_vec;
   std::vector<Eigen::Vector3d> pts/*, dual_pts, face_bars*/;
   std::vector<std::vector<uint32_t>> E_cluster,P_cluster;
   
   int32_t sgn(int32_t val)
   {
      return (0 < val) - (val < 0);
   }
   
   Eigen::Vector3d GetElectricField(uint32_t cube)
   {
	  std::vector<T> u;
	  
	  for (uint32_t i=0; i<12; i++)
		  u.push_back(U[E_cluster[cube][i]]);
	  
      Eigen::Vector3d ret = (u[0]*dual_area_x +	u[ 1]*dual_area_y + u[ 2]*dual_area_z +
	                         u[3]*dual_area_y + u[ 4]*dual_area_z + u[ 5]*dual_area_x +
	                         u[6]*dual_area_z +	u[ 7]*dual_area_z + u[ 8]*dual_area_x +
	                         u[9]*dual_area_y +	u[10]*dual_area_y + u[11]*dual_area_x)/volume;
      return ret;
   }

   Eigen::Vector3d GetMagneticField(uint32_t cube)
   {
	  std::vector<T> u;
	  
	  for (uint32_t i=0; i<6; i++)
		  u.push_back(F[D[cube][i]]);
	  
      Eigen::Vector3d ret = (u[0]*area_z_vec + u[1]*area_y_vec + u[2]*area_x_vec +
	                         u[3]*area_x_vec + u[4]*area_y_vec + u[5]*area_z_vec)
							 /mu[material[cube]]/volume;
      return ret;
   }
   
   Eigen::Vector3d face_barycenter(uint32_t index)
   {
	  std::set<uint32_t> nodes;
	  std::vector<Eigen::Vector3d> ret;
	  for (auto e : C[index])
         for (auto n : G[e])
			 nodes.insert(n);
      for (auto n : nodes)
         ret.push_back(pts[n]);
      return VectorAverage(ret);
   }
   
   Eigen::Vector3d edge_barycenter(uint32_t index)
   {
      return 0.5*(pts[G[index][0]]+pts[G[index][1]]);
   }
   
   Eigen::Vector3d edge_vector(uint32_t index)
   {
      return pts[G[index][1]]-pts[G[index][0]];
   } 
   
   Eigen::Vector3d VectorAverage(std::vector<Eigen::Vector3d> vec)
   {
      Eigen::Vector3d ret(0,0,0);
	  for (auto v : vec)
		  ret += v;
	  return (1/double(vec.size()))*ret;
   }
   
   bool IsGridVol(T x, T y, T z)
   {
	   if (IsGridPoint(x,y,z))
		   return true;
	   if (IsGridPoint(x+Lx,y,z))
		   return true;
	   if (IsGridPoint(x,y+Ly,z))
		   return true;
	   if (IsGridPoint(x,y,z+Lz))
		   return true;
	   if (IsGridPoint(x+Lx,y+Ly,z))
		   return true;
	   if (IsGridPoint(x+Lx,y,z+Lz))
		   return true;
	   if (IsGridPoint(x,y+Ly,z+L))
		   return true;
	   if (IsGridPoint(x+Lx,y+Ly,z+Lz))
		  return true;
	   return false;
   }
   
   bool IsGridPoint(T x, T y, T z)
   {
      //check if it satisfies the set of inequalities
      if (!( x>=xmin ))
         return false;
      if (!( x<=xmax ))
         return false;
      if (!( y>=ymin ))
         return false;
      if (!( y<=ymax ))
         return false;
      if (!( z>=zmin ))
         return false;
      if (!( z<=zmax ))
         return false;
      return true;
   }
};
#endif
