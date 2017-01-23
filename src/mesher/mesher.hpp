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
#include <thread>
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
	   int32_t nv,nf,ne,np;
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
      L=0.0005555555555555555555555555555555555556;
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
	  
	  Eigen::Vector3d dual_area_z(0,0,Lx*Ly/4);
	  Eigen::Vector3d dual_area_y(0,Lx*Lz/4,0);
	  Eigen::Vector3d dual_area_x(Ly*Lz/4,0,0);
	  
	  this->dual_area_x=dual_area_x;
	  this->dual_area_y=dual_area_y;
	  this->dual_area_z=dual_area_z;
	  
	  uint32_t Nx = (fabs(xmax-xmin)) / Lx;// + 1;
	  uint32_t Ny = (fabs(ymax-ymin)) / Ly;// + 1;
	  uint32_t Nz = (fabs(zmax-zmin)) / Lz;// + 1;
	  
	  tot_E= Nx*(Ny+1)*(Nz+1)+(Nx+1)*Ny*(Nz+1)+(Nx+1)*(Ny+1)*Nz;
	  tot_F= Nx*Ny*(Nz+1)+Nz*Nx*(Ny+1)+Ny*Nz*(Nx+1);
	  
	  typedef Eigen::Triplet<uint32_t> Q;
      
	  Eigen::SparseMatrix<uint32_t> previous_layer(Nx,Ny);
	  std::vector<T> average_eps(tot_E,0), average_ni(tot_F,0), face_area(tot_F), edge_len(tot_E);
	  std::vector<uint8_t> boundary_face(tot_F,1); 
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
						std::cout << "# of element to be probed: " << nv+1 << std::endl;
					 }
                  }
			      nv++;
				  std::vector<int32_t> dummy(6), dummyf;
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
                        D.push_back(std::vector<int32_t>({-(nf+1),-(nf+2),-(nf+3),
						                                      nf+4,nf+5,nf+6}));
						nf+=6;
						
						for (uint32_t cnt=0; cnt<6; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						
						E_cluster.push_back(std::vector<int32_t>({ne+1,ne+2,ne+3,ne+4,ne+5,ne+6,
						                                    ne+7,ne+8,ne+9,ne+10,ne+11,ne+12}));
						ne+=12;
						
						for (uint32_t cnt=0; cnt<12; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						
						P_cluster.push_back(std::vector<int32_t>({np+1,np+2,np+3,np+4,
						                                    np+5,np+6,np+7,np+8}));
						np+=8;
						
						for (uint32_t cnt=0; cnt<8; cnt++)
                           Gt.push_back(dummyf);
						
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
						D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(nf+1),-(nf+2),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],ne+1,E_cluster[bottom-1][10],
						                                    ne+2,E_cluster[bottom-1][11],ne+3,ne+4,ne+5,ne+6,ne+7,ne+8}));
						ne+=8;
						
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}

						P_cluster.push_back(std::vector<int32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    np+1,np+2,np+3,np+4}));
						np+=4;
 						for (uint32_t cnt=0; cnt<4; cnt++)
							Gt.push_back(dummyf);
						
                        pts.push_back(pp+inc_z);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 2 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<int32_t>({-(nf+1),-(D[left-1][4]),-(nf+2),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({E_cluster[left-1][5],ne+1,E_cluster[left-1][6],ne+2,E_cluster[left-1][7],ne+3,
						                                    ne+4,ne+5,E_cluster[left-1][11],ne+6,ne+7,ne+8}));
						ne+=8;
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<int32_t>({P_cluster[left-1][2],P_cluster[left-1][3],np+1,np+2,
						                                    P_cluster[left-1][6],P_cluster[left-1][7],np+3,np+4}));
						np+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
							Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_y);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 3 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
                        D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(D[left-1][4]),-(nf+1),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],E_cluster[left-1][6],E_cluster[bottom-1][10],E_cluster[left-1][7],E_cluster[bottom-1][11],
						                                    ne+1,ne+2,E_cluster[left-1][11],ne+3,ne+4,ne+5}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<int32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    P_cluster[left-1][6],P_cluster[left-1][7],np+1,np+2}));
						np+=2;
						for (uint32_t cnt=0; cnt<2; cnt++)
							Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 4 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<int32_t>({-(nf+1),-(nf+2),-(D[back-1][3]),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({ne+1,E_cluster[back-1][3],E_cluster[back-1][4],ne+2,ne+3,ne+4,
						                                    E_cluster[back-1][7],ne+5,ne+6,E_cluster[back-1][10],ne+7,ne+8}));
						ne+=8;
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<int32_t>({P_cluster[back-1][1],np+1,P_cluster[back-1][3],np+2,
						                                    P_cluster[back-1][5],np+3,P_cluster[back-1][7],np+4}));
						np+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
							Gt.push_back(dummyf);
						pts.push_back(pp+inc_x);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 5 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(nf+1),-(D[back-1][3]),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],E_cluster[back-1][4],E_cluster[bottom-1][10],ne+1,E_cluster[bottom-1][11],
						                                    E_cluster[back-1][7],ne+2,ne+3,E_cluster[back-1][10],ne+4,ne+5}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<int32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    P_cluster[back-1][5],np+1,P_cluster[back-1][7],np+2}));
						np+=2;
						for (uint32_t cnt=0; cnt<2; cnt++)
							Gt.push_back(dummyf);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 6 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
						D.push_back(std::vector<int32_t>({-(nf+1),-(D[left-1][4]),-(D[back-1][3]),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						for (uint32_t cnt=0; cnt<4; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({E_cluster[left-1][5],E_cluster[back-1][3],E_cluster[left-1][6],ne+1,E_cluster[left-1][7],ne+2,
						                                    E_cluster[back-1][7],ne+3,E_cluster[left-1][11],E_cluster[back-1][10],ne+4,ne+5}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<int32_t>({P_cluster[left-1][2],P_cluster[left-1][3],P_cluster[back-1][3],np+1,
						                                    P_cluster[left-1][6],P_cluster[left-1][7],P_cluster[back-1][7],np+2}));
						np+=2;
						for (uint32_t cnt=0; cnt<2; cnt++)
							Gt.push_back(dummyf);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
					 case 7 :
					 {
                        // std::cout << nv << " --> " << i << " " << j << " " << k << " --> " << bottom << " " << left << " " << back << std::endl;
                        D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(D[left-1][4]),-(D[back-1][3]),
						                                      nf+1,nf+2,nf+3}));
						nf+=3;
						for (uint32_t cnt=0; cnt<3; cnt++)
						{
							Dt.push_back(dummyf);
							C.push_back(dummyf);
							curl.push_back(1);
							face_bars.push_back(dummy_vec);
						}
						E_cluster.push_back(std::vector<int32_t>({E_cluster[bottom-1][8],E_cluster[bottom-1][9],E_cluster[left-1][6],E_cluster[bottom-1][10],E_cluster[left-1][7],E_cluster[bottom-1][11],
						                                    E_cluster[back-1][7],ne+1,E_cluster[left-1][11],E_cluster[back-1][10],ne+2,ne+3}));
						ne+=3;
						for (uint32_t cnt=0; cnt<3; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						P_cluster.push_back(std::vector<int32_t>({P_cluster[bottom-1][4],P_cluster[bottom-1][5],P_cluster[bottom-1][6],P_cluster[bottom-1][7],
						                                    P_cluster[left-1][6],P_cluster[left-1][7],P_cluster[back-1][7],np+1}));
						np+=1;
						Gt.push_back(dummyf);
						
						pts.push_back(pp+inc_z+inc_x+inc_y);
						break;
					 }
				  }
				  // std::cout << Dt.size() << " --> " << abs(D[nv-1][0]) << " " << abs(D[nv-1][1]) << " " << abs(D[nv-1][2])
				            // << " " << abs(D[nv-1][3]) << " " << abs(D[nv-1][4]) << " " << abs(D[nv-1][5]) << std::endl;
                  Dt[abs(D[nv-1][0])-1].push_back(-nv); 
                  Dt[abs(D[nv-1][1])-1].push_back(-nv);
				  Dt[abs(D[nv-1][2])-1].push_back(-nv);
				  Dt[abs(D[nv-1][3])-1].push_back( nv);
				  Dt[abs(D[nv-1][4])-1].push_back( nv);
				  Dt[abs(D[nv-1][5])-1].push_back( nv);

				  material.push_back(1);
				
				 if (!G[E_cluster[nv-1][0]-1].size())
				 {
					G[E_cluster[nv-1][0]-1]  = std::vector<int32_t>({-P_cluster[nv-1][0],P_cluster[nv-1][1]});
					edge_len[E_cluster[nv-1][0]-1] = Lx;
					Gt[P_cluster[nv-1][0]-1].push_back(-E_cluster[nv-1][0]);
					Gt[P_cluster[nv-1][1]-1].push_back( E_cluster[nv-1][0]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][1]-1].size())
				 {
					G[E_cluster[nv-1][1]-1]  = std::vector<int32_t>({-P_cluster[nv-1][0],P_cluster[nv-1][2]});
					edge_len[E_cluster[nv-1][1]-1] = Ly;
					Gt[P_cluster[nv-1][0]-1].push_back(-E_cluster[nv-1][1]);
					Gt[P_cluster[nv-1][2]-1].push_back( E_cluster[nv-1][1]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][2]-1].size())
				 {
					G[E_cluster[nv-1][2]-1]  = std::vector<int32_t>({-P_cluster[nv-1][0],P_cluster[nv-1][4]});
					edge_len[E_cluster[nv-1][2]-1] = Lz;
					Gt[P_cluster[nv-1][0]-1].push_back(-E_cluster[nv-1][2]);
					Gt[P_cluster[nv-1][4]-1].push_back( E_cluster[nv-1][2]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][3]-1].size())
				 {
					G[E_cluster[nv-1][3]-1]  = std::vector<int32_t>({-P_cluster[nv-1][1],P_cluster[nv-1][3]});
					edge_len[E_cluster[nv-1][3]-1] = Ly;
					Gt[P_cluster[nv-1][1]-1].push_back(-E_cluster[nv-1][3]);
					Gt[P_cluster[nv-1][3]-1].push_back( E_cluster[nv-1][3]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][4]-1].size())
				 {
					G[E_cluster[nv-1][4]-1]  = std::vector<int32_t>({-P_cluster[nv-1][1],P_cluster[nv-1][5]});
					edge_len[E_cluster[nv-1][4]-1] = Lz;
					Gt[P_cluster[nv-1][1]-1].push_back(-E_cluster[nv-1][4]);
					Gt[P_cluster[nv-1][5]-1].push_back( E_cluster[nv-1][4]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][5]-1].size())
				 {
					G[E_cluster[nv-1][5]-1]  = std::vector<int32_t>({-P_cluster[nv-1][2],P_cluster[nv-1][3]});
					edge_len[E_cluster[nv-1][5]-1] = Lx;
					Gt[P_cluster[nv-1][2]-1].push_back(-E_cluster[nv-1][5]);
					Gt[P_cluster[nv-1][3]-1].push_back( E_cluster[nv-1][5]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][6]-1].size())
				 {
					G[E_cluster[nv-1][6]-1]  = std::vector<int32_t>({-P_cluster[nv-1][2],P_cluster[nv-1][6]});
					edge_len[E_cluster[nv-1][6]-1] = Lz;
					Gt[P_cluster[nv-1][2]-1].push_back(-E_cluster[nv-1][6]);
					Gt[P_cluster[nv-1][6]-1].push_back( E_cluster[nv-1][6]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][7]-1].size())
				 {
					G[E_cluster[nv-1][7]-1]  = std::vector<int32_t>({-P_cluster[nv-1][3],P_cluster[nv-1][7]});
					edge_len[E_cluster[nv-1][7]-1] = Lz;
					Gt[P_cluster[nv-1][3]-1].push_back(-E_cluster[nv-1][7]);
					Gt[P_cluster[nv-1][7]-1].push_back( E_cluster[nv-1][7]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][8]-1].size())
				 {
					G[E_cluster[nv-1][8]-1]  = std::vector<int32_t>({-P_cluster[nv-1][4],P_cluster[nv-1][5]});
					edge_len[E_cluster[nv-1][8]-1] = Lx;
					Gt[P_cluster[nv-1][4]-1].push_back(-E_cluster[nv-1][8]);
					Gt[P_cluster[nv-1][5]-1].push_back( E_cluster[nv-1][8]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][9]-1].size())
				 {
					G[E_cluster[nv-1][9]-1]  = std::vector<int32_t>({-P_cluster[nv-1][4],P_cluster[nv-1][6]});
					edge_len[E_cluster[nv-1][9]-1] = Ly;
					Gt[P_cluster[nv-1][4]-1].push_back(-E_cluster[nv-1][9]);
					Gt[P_cluster[nv-1][6]-1].push_back( E_cluster[nv-1][9]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][10]-1].size())
				 {
					G[E_cluster[nv-1][10]-1]  = std::vector<int32_t>({-P_cluster[nv-1][5],P_cluster[nv-1][7]});
					edge_len[E_cluster[nv-1][10]-1] = Ly;
					Gt[P_cluster[nv-1][5]-1].push_back(-E_cluster[nv-1][10]);
					Gt[P_cluster[nv-1][7]-1].push_back( E_cluster[nv-1][10]);
					
					dual_curl.push_back(-1);
					U.push_back(0);
				 }
				 if (!G[E_cluster[nv-1][11]-1].size())
				 {
					G[E_cluster[nv-1][11]-1]  = std::vector<int32_t>({-P_cluster[nv-1][6],P_cluster[nv-1][7]});
					edge_len[E_cluster[nv-1][11]-1] = Lx;
					Gt[P_cluster[nv-1][6]-1].push_back(-E_cluster[nv-1][11]);
					Gt[P_cluster[nv-1][7]-1].push_back( E_cluster[nv-1][11]);
					
					dual_curl.push_back(1);
					U.push_back(0);
				 }
				 // });		 
				 // std::thread C_thread([&] {
				 if (!C[abs(D[nv-1][0])-1].size())
				 {
					C[abs(D[nv-1][0])-1] = std::vector<int32_t>({E_cluster[nv-1][0],-E_cluster[nv-1][1],E_cluster[nv-1][3],-E_cluster[nv-1][5]});
					
					face_area[abs(D[nv-1][0])-1] = area_z;
					// face_bars[abs(D[nv-1][0])-1] = face_barycenter(abs(D[nv-1][0])-1);
					Ct[E_cluster[nv-1][0]-1].push_back( abs(D[nv-1][0]));
					Ct[E_cluster[nv-1][1]-1].push_back(-abs(D[nv-1][0]));
					Ct[E_cluster[nv-1][3]-1].push_back( abs(D[nv-1][0]));
					Ct[E_cluster[nv-1][5]-1].push_back(-abs(D[nv-1][0]));

					F.push_back(0);
				 }
				 else
					 boundary_face[abs(D[nv-1][0])-1]=0;
					
				 
				 if (!C[abs(D[nv-1][1])-1].size())
				 {
					C[abs(D[nv-1][1])-1] = std::vector<int32_t>({-E_cluster[nv-1][0],E_cluster[nv-1][2],-E_cluster[nv-1][4],E_cluster[nv-1][8]});
					curl[abs(D[nv-1][1])-1] = -1;
					face_area[abs(D[nv-1][1])-1] = area_y;
					// face_bars[abs(D[nv-1][1])-1] = face_barycenter(abs(D[nv-1][1])-1);
					Ct[E_cluster[nv-1][0]-1].push_back(-abs(D[nv-1][1]));
					Ct[E_cluster[nv-1][2]-1].push_back( abs(D[nv-1][1]));
					Ct[E_cluster[nv-1][4]-1].push_back(-abs(D[nv-1][1]));
					Ct[E_cluster[nv-1][8]-1].push_back( abs(D[nv-1][1]));

					F.push_back(0);
				 }
				 else
					 boundary_face[abs(D[nv-1][1])-1]=0;
				 
				 if (!C[abs(D[nv-1][2])-1].size())
				 {
					C[abs(D[nv-1][2])-1] = std::vector<int32_t>({E_cluster[nv-1][1],-E_cluster[nv-1][2],E_cluster[nv-1][6],-E_cluster[nv-1][9]});

					face_area[abs(D[nv-1][2])-1] = area_x;
					// face_bars[abs(D[nv-1][2])-1] = face_barycenter(abs(D[nv-1][2])-1);
					Ct[E_cluster[nv-1][1]-1].push_back( abs(D[nv-1][2]));
					Ct[E_cluster[nv-1][2]-1].push_back(-abs(D[nv-1][2]));
					Ct[E_cluster[nv-1][6]-1].push_back( abs(D[nv-1][2]));
					Ct[E_cluster[nv-1][9]-1].push_back(-abs(D[nv-1][2]));
					
					F.push_back(0);
				 }
				 else
					 boundary_face[abs(D[nv-1][2])-1]=0;
				 
				 if (!C[abs(D[nv-1][3])-1].size())
				 {
					C[abs(D[nv-1][3])-1] = std::vector<int32_t>({E_cluster[nv-1][3],-E_cluster[nv-1][4],E_cluster[nv-1][7],-E_cluster[nv-1][10]});

					face_area[abs(D[nv-1][3])-1] = area_x;
					// face_bars[abs(D[nv-1][3])-1] = face_barycenter(abs(D[nv-1][3])-1);
					Ct[E_cluster[nv-1][3]-1].push_back( abs(D[nv-1][3]));
					Ct[E_cluster[nv-1][4]-1].push_back(-abs(D[nv-1][3]));
					Ct[E_cluster[nv-1][7]-1].push_back( abs(D[nv-1][3]));
					Ct[E_cluster[nv-1][10]-1].push_back(-abs(D[nv-1][3]));
					
					F.push_back(0);
				 }
				 else
					 boundary_face[abs(D[nv-1][3])-1]=0;
				 
				 if (!C[abs(D[nv-1][4])-1].size())
				 {
					C[abs(D[nv-1][4])-1] = std::vector<int32_t>({-E_cluster[nv-1][5],E_cluster[nv-1][6],-E_cluster[nv-1][7],E_cluster[nv-1][11]});
					curl[abs(D[nv-1][4])-1] = -1;
					face_area[abs(D[nv-1][4])-1] = area_y;
					// face_bars[abs(D[nv-1][4])-1] = face_barycenter(abs(D[nv-1][4])-1);
					Ct[E_cluster[nv-1][5]-1].push_back(-abs(D[nv-1][4]));
					Ct[E_cluster[nv-1][6]-1].push_back( abs(D[nv-1][4]));
					Ct[E_cluster[nv-1][7]-1].push_back(-abs(D[nv-1][4]));
					Ct[E_cluster[nv-1][11]-1].push_back(abs(D[nv-1][4]));
					
					F.push_back(0);
				 }
				 else
					 boundary_face[abs(D[nv-1][4])-1]=0;
				 
				 if (!C[abs(D[nv-1][5])-1].size())
				 {
					C[abs(D[nv-1][5])-1] = std::vector<int32_t>({E_cluster[nv-1][8],-E_cluster[nv-1][9],E_cluster[nv-1][10],-E_cluster[nv-1][11]});

					face_area[abs(D[nv-1][5])-1] = area_z;
					// face_bars[abs(D[nv-1][5])-1] = face_barycenter(abs(D[nv-1][5])-1);
					Ct[E_cluster[nv-1][8]-1].push_back( abs(D[nv-1][5])); 
					Ct[E_cluster[nv-1][9]-1].push_back(-abs(D[nv-1][5]));
					Ct[E_cluster[nv-1][10]-1].push_back( abs(D[nv-1][5]));
					Ct[E_cluster[nv-1][11]-1].push_back(-abs(D[nv-1][5]));
					
					F.push_back(0);
				 }
				 else
					 boundary_face[abs(D[nv-1][5])-1]=0;
				 
				 
					average_ni[abs(D[nv-1][0])-1] += Lz/2/mu[material[nv-1]]/face_area[abs(D[nv-1][0])-1];
					average_ni[abs(D[nv-1][1])-1] += Ly/2/mu[material[nv-1]]/face_area[abs(D[nv-1][1])-1];
					average_ni[abs(D[nv-1][2])-1] += Lx/2/mu[material[nv-1]]/face_area[abs(D[nv-1][2])-1];
					average_ni[abs(D[nv-1][3])-1] += Lx/2/mu[material[nv-1]]/face_area[abs(D[nv-1][3])-1];
					average_ni[abs(D[nv-1][4])-1] += Ly/2/mu[material[nv-1]]/face_area[abs(D[nv-1][4])-1];
					average_ni[abs(D[nv-1][5])-1] += Lz/2/mu[material[nv-1]]/face_area[abs(D[nv-1][5])-1];
					
					average_eps[E_cluster[nv-1][ 0]-1] += da_x*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 0]-1];
					average_eps[E_cluster[nv-1][ 1]-1] += da_y*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 1]-1];
					average_eps[E_cluster[nv-1][ 2]-1] += da_z*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 2]-1];
					average_eps[E_cluster[nv-1][ 3]-1] += da_y*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 3]-1];
					average_eps[E_cluster[nv-1][ 4]-1] += da_z*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 4]-1];
					average_eps[E_cluster[nv-1][ 5]-1] += da_x*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 5]-1];
					average_eps[E_cluster[nv-1][ 6]-1] += da_z*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 6]-1];
					average_eps[E_cluster[nv-1][ 7]-1] += da_z*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 7]-1];
					average_eps[E_cluster[nv-1][ 8]-1] += da_x*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 8]-1];
					average_eps[E_cluster[nv-1][ 9]-1] += da_y*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][ 9]-1];
					average_eps[E_cluster[nv-1][10]-1] += da_y*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][10]-1];
					average_eps[E_cluster[nv-1][11]-1] += da_x*epsilon[material[nv-1]]/edge_len[E_cluster[nv-1][11]-1];
					
					dual_pts.push_back(pp+0.5*inc_x+0.5*inc_y+0.5*inc_z);
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
	  
      t_mesh.toc();
	  std::cout << "elements: " << nv << " surfaces: " << nf << " edges: " << ne << " nodes: " << np << " dofs: " << nf+ne << std::endl;
	  std::cout << "Meshing takes " << t_mesh << " seconds!" << std::endl;
	  // std::cout << pts[pts.size()-1] << std::endl;
	  
	  t_mesh.tic();

	  M_ni = average_ni;
	  
	  for (uint32_t i=0; i<ne; i++)
	  {
         M_h.push_back(1/average_eps[i]);
		 uint8_t in_b=0;
		 for (auto ff : Ct[i])
		 {
			auto abs_ff = abs(ff)-1;
			if (boundary_face[abs_ff])
			{
				in_b = 1;
				this->is_boundary.push_back(true);
				break;
			}
		 }
		 
		 switch (in_b)
		 {
			case  1 :
			{
				// this->is_boundary.push_back(true);

				auto bar = edge_barycenter(i);
				if (bar(2) <= zmin)
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
	  std::cout << "Constructing material relations takes " << t_mesh << " seconds!" << std::endl;
   }
   
   bool Run(double simulation_time)
   {
	  double step_time_average=0;
	  uint32_t N_of_steps=simulation_time/t_step;
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
			if (is_boundary[j])
				U[j]=time_function*bc[j];
			else
				U[j] = U[j] + t_step*M_h[j]*dual_curl[j]*(F[abs(Ct[j][0])-1]-F[abs(Ct[j][1])-1]+F[abs(Ct[j][2])-1]-F[abs(Ct[j][3])-1]);
		 }

         for (size_t j=0; j<F.size(); j++)
			F[j] = F[j] - t_step*M_ni[j]*curl[j]*(U[abs(C[j][0])-1]-U[abs(C[j][1])-1]+U[abs(C[j][2])-1]-U[abs(C[j][3])-1]);
		
		 auto num_val = GetElectricfield(probe_elem);
		 numeric_values.push_back(num_val(1));
		 numeric_times.push_back(i*t_step);
		 step_cost.toc();
		 step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
		 
		 if ((i+1) % 140 == 0)
			std::cout << "Progress: " << 100*i/N_of_steps << "% done in " << step_time_average << "s, " << step_time_average/i << " s/step" << std::endl;
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
   
   private:
   T xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz,L, volume;
   uint32_t probe_elem, tot_E, tot_F;
   std::vector<T> bc;
   std::vector<bool> is_boundary;
   std::vector<T> M_ni, M_h, U, F;
   std::map<uint32_t,double> epsilon,mu;
   double t_step, freq;
   std::vector<uint32_t> material;
   std::vector<std::vector<int32_t>> D,Dt,C,Ct,G,Gt;
   std::vector<int32_t> curl, dual_curl;
   Eigen::Vector3d dual_area_z, dual_area_y, dual_area_x;
   std::vector<Eigen::Vector3d> pts, dual_pts, face_bars;
   std::vector<std::vector<int32_t>> E_cluster,P_cluster;
   
   int32_t sgn(int32_t val)
   {
      return (0 < val) - (val < 0);
   }
   
   Eigen::Vector3d GetElectricfield(uint32_t cube)
   {
	  std::vector<T> u;
	  
	  for (uint32_t i=0; i<12; i++)
		  u.push_back(U[E_cluster[cube][i]-1]);
	  
      Eigen::Vector3d ret = (u[0]*dual_area_x +	u[ 1]*dual_area_y + u[ 2]*dual_area_z +
	                         u[3]*dual_area_y + u[ 4]*dual_area_z + u[ 5]*dual_area_x +
	                         u[6]*dual_area_z +	u[ 7]*dual_area_z + u[ 8]*dual_area_x +
	                         u[9]*dual_area_y +	u[10]*dual_area_y + u[11]*dual_area_x)/volume;
      return ret;
   }
   
   Eigen::Vector3d face_barycenter(uint32_t index)
   {
	  std::set<uint32_t> nodes;
	  std::vector<Eigen::Vector3d> ret;
	  for (auto e : C[index])
         for (auto n : G[abs(e)-1])
			 nodes.insert(abs(n)-1);
      for (auto n : nodes)
         ret.push_back(pts[n]);
      return VectorAverage(ret);
   }
   
   Eigen::Vector3d edge_barycenter(uint32_t index)
   {
      return 0.5*(pts[abs(G[index][0])-1]+pts[abs(G[index][1])-1]);
   }
   
   Eigen::Vector3d edge_vector(uint32_t index)
   {
      return pts[abs(G[index][1])-1]-pts[abs(G[index][0])-1];
   } 
   
   Eigen::Vector3d VectorAverage(std::vector<Eigen::Vector3d> vec)
   {
      Eigen::Vector3d ret(0,0,0);
	  for (auto v : vec)
		  ret += v;
	  return (1/double(vec.size()))*ret;
   }
   
   // double edge_len(uint32_t index)
   // {
      // return (pts[abs(G[index][0])-1]-pts[abs(G[index][1])-1]).norm();
   // }
   
   // double face_area(uint32_t index)
   // {
      // return edge_len(abs(C[index][0])-1)*edge_len(abs(C[index][1])-1);
   // }
   
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
