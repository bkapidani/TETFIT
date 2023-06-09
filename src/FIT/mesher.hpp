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
#include "Utilities.hpp"
#include <silo.h>
#include <string.h>
#include <set>
#include <cassert>
#include <fstream>
#include <iomanip>

/*general parameters*/
// double pi = 3.141592653589793;
// double mu0 = 4*pi*1e-7;
// double epsilon0 = 8.854187817e-12;
// double c0 = 1 / sqrt( mu0 * epsilon0 );


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
	   
      //Numerical limits
		 
      xmin= 0;
      xmax= 0.05;
      ymin= 0;
      ymax= 0.025;
      zmin= 0;
      zmax= 0.1;
      L=0.0005;
      epsilon[1]=epsilon0;
      mu[1]=mu0;
      sigma[1]=0;
      mag_sigma[1]=0;
      epsilon[2]=epsilon0;
      mu[2]=mu0;
      sigma[2]=0;
      mag_sigma[2]=0;
      freq=5e9;
      xstart = 0.0001;
      xstep  = 0.0001;
      xstop  = 0.0499;
      ystart = 0.0125;
      ystep  = 1;
      ystop  = 0.0125;
      zstart = 0.05;
      zstep  = 1;
      zstop  = 0.05;
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
	  
	assert(xstart <= xstop);
	assert(xstep > 0);
	assert(ystart <= ystop);
	assert(ystep > 0); 
	assert(zstart <= zstop);
	assert(zstep > 0);
	
	 std::vector<Eigen::Vector3d> probe_points;
	double xiter, yiter, ziter;
	for (ziter=zstart; ziter<=zstop; ziter+=zstep)
		for (yiter=ystart; yiter<=ystop; yiter+=ystep)
			for (xiter=xstart; xiter<=xstop; xiter+=xstep)
				probe_points.push_back(Eigen::Vector3d({xiter,yiter,ziter}));
	  
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
	  std::vector<T> average_eps(tot_E,0), average_sigma(tot_E,0), edge_len(tot_E); 
	  std::vector<T> average_ni(tot_F,0), average_mag_sigma(tot_F,0), face_area(tot_F);
	  std::vector<uint8_t> is_ele_lossy(tot_E,0), is_mag_lossy(tot_F,0);
	  // previous_layer.SetZero();
	
      Eigen::Vector3d inc_x(Lx,0,0), inc_y(0,Ly,0), inc_z(0,0,Lz), dummy_vec;
	  std::vector<bool> not_found(probepoints.size(),true);
	  std::vector<int32_t> cplus=std::vector<int32_t>({1,-1,1,-1});
	  std::vector<int32_t> cminus=std::vector<int32_t>({-1,1,-1,1});
	  
	  for(uint32_t k=0;k<Nz;k++)
      {
         py=ymin;
		 std::vector<uint32_t> old_col;
		 Eigen::SparseMatrix<uint32_t> this_layer(Nx,Ny);
         std::vector<Q> tripletList;
         tripletList.reserve(Nx*Ny);

         for (uint32_t j=0;j<Ny;++j)
		 {
            std::vector<uint32_t> this_col(Nx,0);
            px=xmin;
            for (uint32_t i=0;i<Nx;++i)
			{
			   
               if (IsGridVol(px,py,pz))
				   material.push_back(2);
			   else
				   material.push_back(1);
				  Eigen::Vector3d pp(px,py,pz);
				  auto pp2 = pp+inc_x+inc_y+inc_z;
				  
				  for (uint32_t piter=0; piter<probe_points.size(); ++piter)
				  {
					  if (not_found[piter])
					  {
						 
						 if (probe_points[piter](0)>= pp(0) && probe_points[piter](0)<= pp2(0) && 
							 probe_points[piter](1)>= pp(1) && probe_points[piter](1)<= pp2(1) &&
							 probe_points[piter](2)>= pp(2) && probe_points[piter](2)<= pp2(2))
						 {
							not_found[piter] = false;
							probe_elem.push_back(nv);
							probepoints.push_back(probe_points[piter]);
							// std::cout << "# of element to be probed: " << nv+1 << std::endl;
						 }
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

				  // material.push_back(1);
				
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

					average_ni[D[nv][0]] += Lz/2/mu[material[nv]];
					average_ni[D[nv][1]] += Ly/2/mu[material[nv]];
					average_ni[D[nv][2]] += Lx/2/mu[material[nv]];
					average_ni[D[nv][3]] += Lx/2/mu[material[nv]];
					average_ni[D[nv][4]] += Ly/2/mu[material[nv]];
					average_ni[D[nv][5]] += Lz/2/mu[material[nv]];
					
					if (mag_sigma[material[nv]] != 0)
					{
						average_mag_sigma[D[nv][0]] += Lz/2/mag_sigma[material[nv]]; is_mag_lossy[D[nv][0]]++;
						average_mag_sigma[D[nv][1]] += Ly/2/mag_sigma[material[nv]]; is_mag_lossy[D[nv][1]]++;
						average_mag_sigma[D[nv][2]] += Lx/2/mag_sigma[material[nv]]; is_mag_lossy[D[nv][2]]++;
						average_mag_sigma[D[nv][3]] += Lx/2/mag_sigma[material[nv]]; is_mag_lossy[D[nv][3]]++;
						average_mag_sigma[D[nv][4]] += Ly/2/mag_sigma[material[nv]]; is_mag_lossy[D[nv][4]]++;
						average_mag_sigma[D[nv][5]] += Lz/2/mag_sigma[material[nv]]; is_mag_lossy[D[nv][5]]++;
					}
					
					average_eps[E_cluster[nv][ 0]] += da_x*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 1]] += da_y*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 2]] += da_z*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 3]] += da_y*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 4]] += da_z*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 5]] += da_x*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 6]] += da_z*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 7]] += da_z*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 8]] += da_x*epsilon[material[nv]];
					average_eps[E_cluster[nv][ 9]] += da_y*epsilon[material[nv]];
					average_eps[E_cluster[nv][10]] += da_y*epsilon[material[nv]];
					average_eps[E_cluster[nv][11]] += da_x*epsilon[material[nv]];
					
					if (sigma[material[nv]] != 0)
					{
						average_sigma[E_cluster[nv][ 0]] += da_x*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 0]]++;
						average_sigma[E_cluster[nv][ 1]] += da_y*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 1]]++;
						average_sigma[E_cluster[nv][ 2]] += da_z*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 2]]++;
						average_sigma[E_cluster[nv][ 3]] += da_y*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 3]]++;
						average_sigma[E_cluster[nv][ 4]] += da_z*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 4]]++;
						average_sigma[E_cluster[nv][ 5]] += da_x*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 5]]++;
						average_sigma[E_cluster[nv][ 6]] += da_z*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 6]]++;
						average_sigma[E_cluster[nv][ 7]] += da_z*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 7]]++;
						average_sigma[E_cluster[nv][ 8]] += da_x*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 8]]++;
						average_sigma[E_cluster[nv][ 9]] += da_y*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][ 9]]++;
						average_sigma[E_cluster[nv][10]] += da_y*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][10]]++;
						average_sigma[E_cluster[nv][11]] += da_x*sigma[material[nv]]; is_ele_lossy[E_cluster[nv][11]]++;
					}
					
					// dual_pts.push_back(pp+0.5*inc_x+0.5*inc_y+0.5*inc_z);
					nv++;
					tripletList.push_back(Q(i,j,nv));
					this_col[i]=nv;
			   // }
			   
               px+=Lx;
			}
			
			old_col=std::move(this_col);
			py+=Ly;
		 }

         this_layer.setFromTriplets(tripletList.begin(), tripletList.end());
		 previous_layer=std::move(this_layer);
		 pz+=Lz;
      }
	  
			Nvec.resize(tot_F);
			Hvec.resize(tot_E);
			

	  for (uint32_t i=0; i<nf; ++i)
	  {
			Nvec(i)=average_ni[i]/face_area[i];

		if (is_mag_lossy[i])
		{
			M_ni.push_back(1/(face_area[i]*(1/average_ni[i]+0.5*t_step/average_mag_sigma[i])));
			M_mu.push_back(face_area[i]*(1/average_ni[i]-0.5*t_step/average_mag_sigma[i]));
		}
		else
		{
			M_ni.push_back(average_ni[i]/face_area[i]);
			M_mu.push_back(face_area[i]/average_ni[i]);
		}
	  }
	  // uint32_t number_of_lossy=0;
	  
	  for (uint32_t i=0; i<ne; ++i)
	  {
			Hvec(i)=edge_len[i]/average_eps[i];

		 if (is_ele_lossy[i])
		 {
			// std::cout << ++number_of_lossy << std::endl;
			M_h.push_back(edge_len[i]/(average_eps[i] + 0.5*t_step*average_sigma[i]));
			M_e.push_back((average_eps[i] - 0.5*t_step*average_sigma[i])/edge_len[i]);
		 }
		 else
		 {
			M_h.push_back(edge_len[i]/average_eps[i]);
			M_e.push_back(average_eps[i]/edge_len[i]);
		 }
		 
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
					bc.push_back(sin( bar(0)*PI/0.05 )*edge_vector(i)(1));
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
	  
	  this->is_mag_lossy=std::move(is_mag_lossy);
	  this->is_ele_lossy=std::move(is_ele_lossy);
	  
      t_mesh.toc();
	  std::cout << "Meshing and material modeling done in " << t_mesh << " seconds" << std::endl;

   }
   
	/*void RunEigenBased(const double simulation_time)
	{
		std::cout << "---------------------- Running FDTD simulation ----------------------" << std::endl << std::endl;
		double step_time_average=0;
		const uint32_t N_of_steps=simulation_time/t_step;
		uint32_t i,e1,e2,e3,f1,f2,f3;

		// std::cout << Nx << " " << Ny << " " << Nz << std::endl;

		T time_function;
		timecounter step_cost, debug_counter;
		std::vector<T> probe_numeric_yvalues,numeric_times;
		uint32_t nv = D.size()-1;
		typedef Eigen::Triplet<double,uint32_t> W;
		Eigen::SparseMatrix<double> Cmat(C.size(),Ct.size()),CTmat(Ct.size(),C.size());
		std::vector<W> tripletList1, tripletList2;
		tripletList1.reserve(4*C.size());
		tripletList2.reserve(4*Ct.size());
		
		Eigen::VectorXd U_eigen(U.size()), F_eigen(F.size());
		// Eigen::VectorXd Hvec(M_h.data()), Nvec(M_ni.data());
		
		// std::cout << Hvec.size() << "------" << Nvec.size() << std::endl;
		
		for (uint32_t j=0; j<nv; ++j)
		{
			e1 = E_cluster[j][ 7];
			e2 = E_cluster[j][10];
			e3 = E_cluster[j][11];
			f1 = D[j][3];
			f2 = D[j][4];
			f3 = D[j][5];
			
			if (!boundary_face[f1])
			{
				tripletList1.push_back(W(f1,uint32_t(C[f1][0]), double(1)));
				tripletList1.push_back(W(f1,uint32_t(C[f1][1]),-double(1)));
				tripletList1.push_back(W(f1,uint32_t(C[f1][2]), double(1)));
				tripletList1.push_back(W(f1,uint32_t(C[f1][3]),-double(1)));
			}
			
			if (!boundary_face[f2])
			{
				tripletList1.push_back(W(f2,uint32_t(C[f2][0]),-double(1)));
				tripletList1.push_back(W(f2,uint32_t(C[f2][1]), double(1)));
				tripletList1.push_back(W(f2,uint32_t(C[f2][2]),-double(1)));
				tripletList1.push_back(W(f2,uint32_t(C[f2][3]), double(1)));
			}
			
			if (!boundary_face[f3])
			{
				tripletList1.push_back(W(f3,uint32_t(C[f3][0]), double(1)));
				tripletList1.push_back(W(f3,uint32_t(C[f3][1]),-double(1)));
				tripletList1.push_back(W(f3,uint32_t(C[f3][2]), double(1)));
				tripletList1.push_back(W(f3,uint32_t(C[f3][3]),-double(1)));
			}
			
			if (!is_boundary[e1])
			{
				tripletList2.push_back(W(e1,uint32_t(Ct[e1][0]), double(1)));
				tripletList2.push_back(W(e1,uint32_t(Ct[e1][1]),-double(1)));
				tripletList2.push_back(W(e1,uint32_t(Ct[e1][2]), double(1)));
				tripletList2.push_back(W(e1,uint32_t(Ct[e1][3]),-double(1)));
			}
			
			if (!is_boundary[e2])
			{
				tripletList2.push_back(W(e2,uint32_t(Ct[e2][0]),-double(1)));
				tripletList2.push_back(W(e2,uint32_t(Ct[e2][1]), double(1)));
				tripletList2.push_back(W(e2,uint32_t(Ct[e2][2]),-double(1)));
				tripletList2.push_back(W(e2,uint32_t(Ct[e2][3]), double(1)));
			}
			
			if (!is_boundary[e3])
			{
				tripletList2.push_back(W(e3,uint32_t(Ct[e3][0]), double(1)));
				tripletList2.push_back(W(e3,uint32_t(Ct[e3][1]),-double(1)));
				tripletList2.push_back(W(e3,uint32_t(Ct[e3][2]), double(1)));
				tripletList2.push_back(W(e3,uint32_t(Ct[e3][3]),-double(1)));
			}
		}
		
		Cmat.setFromTriplets(tripletList1.begin(), tripletList1.end());
		CTmat.setFromTriplets(tripletList2.begin(), tripletList2.end());
		
		
		for (i=0; i*t_step <= simulation_time; ++i)
		{
			step_cost.tic();
			
			// debug_counter.tic();
			time_function=sin(2*PI*freq*i*t_step);
			for (auto ee : bc_edges)
				if (bc[ee] != 0)
					U_eigen(ee) = time_function*bc[ee];
			
			F_eigen = F_eigen - Nvec.cwiseProduct(Cmat*U_eigen);
			U_eigen = U_eigen + Hvec.cwiseProduct(CTmat*F_eigen);
			// auto num_val = GetElectricField(probe_elem[piter]);
			// probe_numeric_yvalues[piter].push_back(num_val(1));
			// numeric_times.push_back(i*t_step);
			step_cost.toc();
			step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
			
			// debug_counter.toc();
			// std::cout << "Post proc takes " << debug_counter << " seconds" << std::endl;
			// if (i % 20 == 0)
				// ExportFields(i);
			
			if ((i+1) % 140 == 0)
				std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(7) << step_time_average << "s, " 
						  << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
		}
		
		std::ofstream os;
		os.open("numeric_FIT.dat");
		for (uint32_t k=0; k < probe_numeric_yvalues.size(); k++)
		  os << numeric_times[k] << " " << probe_numeric_yvalues[k] << std::endl;
		os.close();
		std::cout << "Time step takes (average) " << step_time_average/(double(i)) << " seconds (" << i << " time steps!)" << std::endl;
		std::cout << "Total running time is "     << step_time_average << " seconds" << std::endl;
	}*/
	
	/*void Run_VoxBased(const double simulation_time)
	{
		std::cout << "---------------------- Running FDTD simulation ----------------------" << std::endl << std::endl;
		double step_time_average=0;
		const uint32_t N_of_steps=simulation_time/t_step;
		uint32_t i,e1,e2,e3,f1,f2,f3;

		// std::cout << Nx << " " << Ny << " " << Nz << std::endl;

		T time_function;
		timecounter step_cost, debug_counter;
		std::vector<std::vector<T>> probe_numeric_yvalues, analytic_values;
		std::vector<T> numeric_times;
		
		uint32_t nv = D.size()-1;
	  
		for (i=0; i*t_step <= simulation_time; ++i)
		{
			step_cost.tic();
			
			// debug_counter.tic();
			time_function=sin(2*PI*freq*i*t_step);
			for (auto ee : bc_edges)
				if (bc[ee] != 0)
					U[ee] = time_function*bc[ee];
			
			// debug_counter.toc();
			// std::cout << "BCs take " << debug_counter << " seconds" << std::endl;
			// debug_counter.tic();
			
			for (uint32_t j=0; j<nv; ++j)
			{
				// if (j==0)
					// debug_counter.tic();
				
				// e1 = E_cluster[j][ 7];
				// e2 = E_cluster[j][10];
				// e3 = E_cluster[j][11];
				// f1 = D[j][3];
				// f2 = D[j][4];
				// f3 = D[j][5];

				// if (j==0)
				// {
					// debug_counter.tic();
				// }
				
				if (!is_boundary[E_cluster[j][ 7]])
				{
					if (is_ele_lossy[E_cluster[j][ 7]]!=0)
						U[E_cluster[j][ 7]] = M_h[E_cluster[j][ 7]]*(M_e[E_cluster[j][ 7]]*U[E_cluster[j][ 7]] + t_step*(F[D[j][3]]-F[D[j+Nx][3]]+F[D[j+1][4]]-F[D[j][4]]));
					else
						U[E_cluster[j][ 7]] += M_h[E_cluster[j][ 7]]*t_step*(F[D[j][3]]-F[D[j+Nx][3]]+F[D[j+1][4]]-F[D[j][4]]);
				}
				if (!is_boundary[E_cluster[j][10]])
				{
					if (is_ele_lossy[E_cluster[j][10]]!=0)
						U[E_cluster[j][10]] = M_h[E_cluster[j][10]]*(M_e[E_cluster[j][10]]*U[E_cluster[j][10]] + t_step*(F[D[j][5]]-F[D[j+1][5]]+F[D[j+Nx*Ny][3]]-F[D[j][3]]));
					else
						U[E_cluster[j][10]] += M_h[E_cluster[j][10]]*t_step*(F[D[j][5]]-F[D[j+1][5]]+F[D[j+Nx*Ny][3]]-F[D[j][3]]);
				}
				if (!is_boundary[E_cluster[j][11]])
				{
					if (is_ele_lossy[E_cluster[j][11]]!=0)
						U[E_cluster[j][11]] = M_h[E_cluster[j][11]]*(M_e[E_cluster[j][11]]*U[E_cluster[j][11]] + t_step*(F[D[j][4]]-F[D[j+Nx*Ny][4]]+F[D[j+Nx][5]]-F[D[j][5]]));
					else
						U[E_cluster[j][11]] += M_h[E_cluster[j][11]]*t_step*(F[D[j][4]]-F[D[j+Nx*Ny][4]]+F[D[j+Nx][5]]-F[D[j][5]]);
				}
				
				// if (j==0)
				// {
					// debug_counter.toc();
					// std::cout << "Electric part takes " << debug_counter << " seconds" << std::endl;
					// debug_counter.tic();
				// }
				
				if (!boundary_face[D[j][3]])
				{
					if (is_mag_lossy[D[j][3]])
						F[D[j][3]] = M_ni[D[j][3]]*(M_mu[D[j][3]]*F[D[j][3]]-t_step*(U[E_cluster[j][ 7]]-U[E_cluster[j][4]]+U[E_cluster[j][3]]-U[E_cluster[j][10]]));
					else
						F[D[j][3]] -= t_step*M_ni[D[j][3]]*(U[E_cluster[j][ 7]]-U[E_cluster[j][4]]+U[E_cluster[j][3]]-U[E_cluster[j][10]]);
				}
				if (!boundary_face[D[j][4]])
				{
					if (is_mag_lossy[D[j][4]])
						F[D[j][4]] = M_ni[D[j][4]]*(M_mu[D[j][4]]*F[D[j][4]]-t_step*(U[E_cluster[j][6]]-U[E_cluster[j][ 7]]+U[E_cluster[j][11]]-U[E_cluster[j][5]]));
					else
						F[D[j][4]] -= t_step*M_ni[D[j][4]]*(U[E_cluster[j][6]]-U[E_cluster[j][ 7]]+U[E_cluster[j][11]]-U[E_cluster[j][5]]);
				}
				if (!boundary_face[D[j][5]])
				{
					if (is_mag_lossy[D[j][5]])
						F[D[j][5]] = M_ni[D[j][5]]*(M_mu[D[j][5]]*F[D[j][5]]-t_step*(U[E_cluster[j][8]]-U[E_cluster[j][11]]+U[E_cluster[j][10]]-U[E_cluster[j][9]]));
					else
						F[D[j][5]] -= t_step*M_ni[D[j][5]]*(U[E_cluster[j][8]]-U[E_cluster[j][11]]+U[E_cluster[j][10]]-U[E_cluster[j][9]]);
				}
				
				// if (j==0)
				// {
					// debug_counter.toc();
					// std::cout << "Magnetic part takes " << debug_counter << " seconds" << std::endl;
					// debug_counter.tic();
				// }
			}
			
			// debug_counter.toc();
			// std::cout << "Other fields take " << debug_counter << " seconds" << std::endl;
			// debug_counter.tic();
			
			auto num_val = GetElectricField(probe_elem[piter]);
			double current_time = i*t_step;
			auto spp = SpaceTimePoint({probepoints[piter](0),probepoints[piter](1),probepoints[piter](2),current_time});
			auto anal_val = analytic_value(spp,sigma[material[probe_elem[piter]]],epsilon[material[probe_elem[piter]]],mu[material[probe_elem[piter]]],freq);
			probe_numeric_yvalues[piter].push_back(num_val(1));
			analytic_values[piter].push_back(anal_val);
			numeric_times.push_back(current_time);
			step_cost.toc();
			step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
			
			// debug_counter.toc();
			// std::cout << "Post proc takes " << debug_counter << " seconds" << std::endl;
			// if (i % 20 == 0)
				// ExportFields(i);
			
			if ((i+1) % 140 == 0)
				std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(7) << step_time_average << "s, " 
						  << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
		}
		
		std::ofstream os;
		os.open("numeric_FIT.dat");
		for (uint32_t k=0; k < probe_numeric_yvalues.size(); k++)
		  os << numeric_times[k] << " " << probe_numeric_yvalues[k] << std::endl;
		os.close();
		std::cout << "Time step takes (average) " << step_time_average/(double(i)) << " seconds (" << i << " time steps!)" << std::endl;
		std::cout << "Total running time is "     << step_time_average << " seconds" << std::endl;
	}*/
   
   void Run(const double simulation_time)
   {
	  std::cout  << std::endl << std::endl <<"---------------- Running FDTD simulation ----------------" << std::endl << std::endl;
	  double step_time_average=0;
	  double max_rel_err = 0;
	  const uint32_t N_of_steps=simulation_time/t_step;
	  uint32_t i;
	  
	  T time_function;
	  timecounter step_cost;
		std::vector<std::vector<T>> probe_numeric_xvalues, probe_numeric_yvalues, probe_numeric_zvalues, analytic_values;
		std::vector<T> numeric_times;
      std::vector<T> U_old(U); std::vector<T> F_old(F);
	  
		std::vector<T> dummy_vector;
			
		for (uint32_t piter=0; piter<probepoints.size(); ++piter)
		{
			probe_numeric_xvalues.push_back(dummy_vector);
			probe_numeric_yvalues.push_back(dummy_vector);
			probe_numeric_zvalues.push_back(dummy_vector);
			analytic_values.push_back(dummy_vector);
		}
		 

		std::cout << std::endl << "Simulation parameters:" 		<< std::endl;
		std::cout << std::setw(20) << "Mesh diamter: " 			<< std::setw(20) << 0.5*L*sqrt(3)				 << "   m" << std::endl;
		std::cout << std::setw(20) << "Simulation time: "  		<< std::setw(20) << simulation_time              << " sec" << std::endl;
		std::cout << std::setw(20) << "Time step: "  			<< std::setw(20) << t_step                       << " sec" << std::endl;
		std::cout << std::setw(20) << "Unknowns: "         		<< std::setw(20) << U.size()+F.size()      		 << std::endl  << std::endl;
      
	  for (i=0; i*t_step <= simulation_time; ++i)
	  {
		 step_cost.tic();
		 time_function=sin(2*PI*freq*(i)*t_step);
		 double current_time = double(i)*t_step;
		 // std::cout << "where does it go" << std::endl;
		 U_old=U;
         for (uint32_t j=0; j<U.size(); ++j)
		 {
			// if (bc[j]!=0)
			// {
				if (!is_boundary[j])
					U[j] = M_h[j]*(M_e[j]*U_old[j] +t_step*dual_curl[j]*(F[Ct[j][0]]-F[Ct[j][1]]+F[Ct[j][2]]-F[Ct[j][3]]));
				else
					U[j] = time_function*bc[j];
			// }
		 }
		 
		 // std::cout << "where does it go" << std::endl;
		 
		 // for (auto j : bc_edges)
			 // U[j] = time_function*bc[j];
		 
		 F_old=F;
         for (uint32_t j=0; j<F.size(); ++j)
			 // if (!boundary_face[j])
				F[j] =  M_ni[j]*(M_mu[j]*F_old[j] - t_step*curl[j]*(U[C[j][0]]-U[C[j][1]]+U[C[j][2]]-U[C[j][3]]));
		
		if ((current_time - t_step) < 1e-9 && (current_time + t_step) > 1e-9)
		{

			
			for (uint32_t piter=0; piter < probepoints.size(); ++piter)
			{
				auto num_val = GetElectricField(probe_elem[piter]);
				double current_time = double(i+1)*t_step;
				auto spp = SpaceTimePoint({probepoints[piter](0),probepoints[piter](1),probepoints[piter](2),current_time});
				auto anal_val = analytic_value(spp,sigma[material[probe_elem[piter]]],epsilon[material[probe_elem[piter]]],mu[material[probe_elem[piter]]],freq);
				// double anal_val = 0;
				probe_numeric_xvalues[piter].push_back(num_val(0));
				probe_numeric_yvalues[piter].push_back(num_val(1));
				probe_numeric_zvalues[piter].push_back(num_val(2));
				analytic_values[piter].push_back(anal_val);

			}
		}
		
		numeric_times.push_back(current_time);
		step_cost.toc();
		step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
		 
		// if ((current_time - t_step) < 1e-9 && (current_time + t_step) > 1e-9)
			// ExportFields(i);
		
		if ((i+1) % 140 == 0)
			std::cout << "-----------" << "Progress: " << 100*i/N_of_steps << "% done in " << std::setw(7) << step_time_average << "s, " 
					  << std::setw(8) << step_time_average/i << std::setw(7) << " s/step" << "-----------" << std::endl;
					  
		// if ((current_time - t_step) < 1e-9 && (current_time + t_step) > 1e-9)
		// {
			// double rel_err = fabs(anal_val - num_val(1))/fabs(anal_val);
			// if (rel_err > max_rel_err)
				// max_rel_err = rel_err;
		// }
	  }
	  
	  /* Output stats and fields*/
	  std::ofstream os, os_a;
	  os.open("numeric_FIT.dat");
	  os_a.open("analytic_FIT.dat");
	  
	  for (uint32_t k=0; k<numeric_times.size(); ++k)
	  {
		  for (uint32_t p=0; p < probepoints.size(); p++)
		  {
					os << std::setw(15) << probepoints[p](0) << " ";
					os << std::setw(15) << probepoints[p](1) << " ";
					os << std::setw(15) << probepoints[p](2) << " ";
					os << std::setw(15) << numeric_times[k] << " ";
					os << std::setw(15) << probe_numeric_xvalues[p][k] << " ";
					os << std::setw(15) << probe_numeric_yvalues[p][k] << " ";
					os << std::setw(15) << probe_numeric_zvalues[p][k];
					os << std::endl;
					
					os_a << std::setw(15) << probepoints[p](0) 			<< " ";
					os_a << std::setw(15) << probepoints[p](1) 			<< " ";
					os_a << std::setw(15) << probepoints[p](2) 			<< " ";
					os_a << std::setw(15) << numeric_times[k] 			<< " ";
					os_a << std::setw(15) <<                          0 << " ";
					os_a << std::setw(15) << analytic_values[p][k] 		<< " ";
					os_a << std::setw(15) <<                           0;
					os_a << std::endl;	  
		  }
	  }
	  
	  os.close(); os_a.close();
	  
	  std::ofstream max_error_file("max_rel_err_fit.dat",std::ofstream::out | std::ofstream::app);
	  max_error_file << U.size()+F.size() << " " << max_rel_err << std::endl;
	  
	  std::cout << "Time step takes (average) " << step_time_average/(double(i)) << " seconds (" << i << " time steps!)" << std::endl;
	  std::cout << "Total running time is "     << step_time_average << " seconds" << std::endl;
	  std::cout << "Max relative error is "     << max_rel_err <<  std::endl;
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
		namedummy << "fit" << time << ".silo";
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
			Hx_vals.push_back(Hfield(0));
			Hy_vals.push_back(Hfield(1));
			Hz_vals.push_back(Hfield(2));
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
		for (uint32_t i=0; i<=Nx; ++i)
			nodex.push_back(i*Lx);
		for (uint32_t i=0; i<=Ny; ++i)
			nodey.push_back(i*Ly);
		for (uint32_t i=0; i<=Nz; ++i)
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
   bool output_to_file = false;
   T xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz,L,volume;
   uint32_t  tot_E, tot_F, Nx, Ny, Nz;
   std::vector<bool> is_boundary;
   std::vector<T> M_ni, M_h, M_e, M_mu, U, F, bc;
   std::map<uint32_t,double> epsilon,mu,sigma,mag_sigma;
   double t_step, freq;
   T xstart, xstep, xstop, ystart, ystep, ystop, zstart, zstep, zstop;
   DBfile *_siloDb=NULL;
   std::vector<uint32_t> probe_elem, material, bc_edges;
   std::vector<uint8_t> boundary_face;
   std::vector<std::vector<uint32_t>> D/*,Dt*/,C,Ct,G/*,Gt*/;
   std::vector<int32_t> curl, dual_curl;
   Eigen::Vector3d dual_area_z, dual_area_y, dual_area_x;
   Eigen::Vector3d area_z_vec, area_y_vec, area_x_vec;
   Eigen::VectorXd Hvec,Nvec;
   std::vector<Eigen::Vector3d> probepoints;
   std::vector<Eigen::Vector3d> pts/*, dual_pts, face_bars*/;
   std::vector<std::vector<uint32_t>> E_cluster,P_cluster;
   std::vector<uint8_t> is_ele_lossy, is_mag_lossy;
   
   int32_t sgn(int32_t val)
   {
      return (0 < val) - (val < 0);
   }
   
   Eigen::Vector3d GetElectricField(uint32_t cube)
   {
	  std::vector<T> u;
	  
	  for (uint32_t i=0; i<12; ++i)
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
	  
	  for (uint32_t i=0; i<6; ++i)
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
	   if (!IsGridPoint(x,y,z))
		   return false;
	   if (!IsGridPoint(x+Lx,y,z))
		   return false;
	   if (!IsGridPoint(x,y+Ly,z))
		   return false;
	   if (!IsGridPoint(x,y,z+Lz))
		   return false;
	   if (!IsGridPoint(x+Lx,y+Ly,z))
		   return false;
	   if (!IsGridPoint(x+Lx,y,z+Lz))
		   return false;
	   if (!IsGridPoint(x,y+Ly,z+L))
		   return false;
	   if (!IsGridPoint(x+Lx,y+Ly,z+Lz))
		  return false;
	   return true;
   }
   
   bool IsGridPoint(T x, T y, T z)
   {
      //check if it satisfies the set of inequalities
      if (!( pow(x-0.025,2) + pow(y-0.0125,2) + pow(z-0.05,2) < pow(0.005,2) ))
         return false;
      return true;
   }
};
#endif
