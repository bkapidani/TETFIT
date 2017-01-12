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

template<typename T>
class mesher
{
   public:
   mesher(void)
   {
	   int32_t nv,nf,ne,np;
	   nv=nf=ne=np=0;
      //Numerical limits
		 
      xmin=-1;
      xmax= 1;
      ymin=-1;
      ymax= 1;
      zmin=-1;
      zmax= 1;
      L=0.1;
      //other instructions in constructor
	  // std::cout << "Qui?" << std::endl;
	  auto px = xmin;
	  auto py = ymin;
	  auto pz = zmin;  
	  
	  Lx = Ly = Lz = L;
      double mu0 = 4*3.141592*1e-7;
      double epsilon0 = 8.854187817e-12;
      double c0 = 1 / sqrt( mu0 * epsilon0 );
	  t_step = 0.5*sqrt(pow(Lx,2)+pow(Ly,2)+pow(Lz,2))/c0;
	  
	  uint32_t Nx = (fabs(xmax-xmin)) / Lx;
	  uint32_t Ny = (fabs(ymax-ymin)) / Ly;
	  uint32_t Nz = (fabs(zmax-zmin)) / Lz;
	  
	  typedef Eigen::Triplet<uint32_t> U;
      std::vector<U> tripletList;
      tripletList.reserve(Nx*Ny);
      
	  Eigen::SparseMatrix<uint32_t> previous_layer(Nx,Ny);
	  
	  // previous_layer.SetZero();
	
      Eigen::Vector3d inc_x(Lx,0,0);
	  Eigen::Vector3d inc_y(0,Ly,0);
	  Eigen::Vector3d inc_z(0,0,Lz);
	  std::vector<std::vector<int32_t>> vte,vtn;
	  
	  for(uint32_t k=0;k<Nz;k++)
      {
         py=ymin;
		 std::vector<uint32_t> old_col;
		 Eigen::SparseMatrix<uint32_t> this_layer(Nx,Ny);
		 
         for (uint32_t j=0;j<Ny;j++)
		 {
            std::vector<uint32_t> this_col(Nx,0);
            px=xmin;
            for (uint32_t i=0;i<Nx;i++)
			{
			   
               if (IsGridVol(px,py,pz))
			   {
				  Eigen::Vector3d pp(px,py,pz);
                  nv++;
				  std::vector<int32_t> dummy(6), dummyf;
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
                        // std::cout << i << " " << j << " " << k << std::endl; 
                        D.push_back(std::vector<int32_t>({-(nf+1),-(nf+2),-(nf+3),
						                                      nf+4,nf+5,nf+6}));
						nf+=6;
						vte.push_back(std::vector<int32_t>({ne+1,ne+2,ne+3,ne+4,ne+5,ne+6,
						                                    ne+7,ne+8,ne+9,ne+10,ne+11,ne+12}));
						ne+=12;
						vtn.push_back(std::vector<int32_t>({np+1,np+2,np+3,np+4,
						                                    np+5,np+6,np+7,np+8}));
						np+=8;
						
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
                        // std::cout << i << " " << j << " " << k << std::endl;
						D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(nf+1),-(nf+2),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],ne+1,vte[bottom-1][10],
						                                    ne+2,vte[bottom-1][11],ne+3,ne+4,ne+5,ne+6,ne+7,ne+8}));
						ne+=8;
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    np+1,np+2,np+3,np+4}));
						np+=4;
 
                        pts.push_back(pp+inc_z);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 2 :
					 {
                        // std::cout << i << " " << j << " " << k << std::endl;
						D.push_back(std::vector<int32_t>({-(nf+1),-(D[left-1][4]),-(nf+2),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						vte.push_back(std::vector<int32_t>({vte[left-1][5],ne+1,vte[left-1][6],ne+2,vte[left-1][7],ne+3,
						                                    ne+4,ne+5,vte[left-1][11],ne+6,ne+7,ne+8}));
						ne+=8;
						vtn.push_back(std::vector<int32_t>({vtn[left-1][2],vtn[left-1][3],np+1,np+2,
						                                    vtn[left-1][6],vtn[left-1][7],np+3,np+4}));
						np+=4;
						
						pts.push_back(pp+inc_y);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 3 :
					 {
                        // std::cout << i << " " << j << " " << k << std::endl;
                        D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(D[left-1][4]),-(nf+1),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[left-1][6],vte[bottom-1][10],vte[left-1][7],vte[bottom-1][11],
						                                    ne+1,ne+2,vte[left-1][11],ne+3,ne+4,ne+5}));
						ne+=5;
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[left-1][6],vtn[left-1][7],np+1,np+2}));
						np+=2;
						 
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 4 :
					 {
                        // std::cout << i << " " << j << " " << k << std::endl;
						D.push_back(std::vector<int32_t>({-(nf+1),-(nf+2),-(D[back-1][3]),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						vte.push_back(std::vector<int32_t>({ne+1,vte[back-1][3],vte[back-1][4],ne+2,ne+3,ne+4,
						                                    vte[back-1][7],ne+5,ne+6,vte[back-1][10],ne+7,ne+8}));
						ne+=8;
						vtn.push_back(std::vector<int32_t>({vtn[back-1][1],np+1,vtn[back-1][3],np+2,
						                                    vtn[back-1][5],np+3,vtn[back-1][7],np+4}));
						np+=4;
						
						pts.push_back(pp+inc_x);
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 5 :
					 {
                        // std::cout << i << " " << j << " " << k << std::endl;
						D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(nf+1),-(D[back-1][3]),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[back-1][4],vte[bottom-1][10],ne+1,vte[bottom-1][11],
						                                    vte[back-1][7],ne+2,ne+3,vte[back-1][10],ne+4,ne+5}));
						ne+=5;
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[back-1][5],np+1,vtn[back-1][7],np+2}));
						np+=2;
						
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 6 :
					 {
                        // std::cout << i << " " << j << " " << k << std::endl;
						D.push_back(std::vector<int32_t>({-(nf+1),-(D[left-1][4]),-(D[back-1][3]),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						vte.push_back(std::vector<int32_t>({vte[left-1][5],vte[back-1][3],vte[left-1][6],ne+1,vte[left-1][7],ne+2,
						                                    vte[back-1][7],ne+3,vte[left-1][11],vte[back-1][10],ne+4,ne+5}));
						ne+=5;
						vtn.push_back(std::vector<int32_t>({vtn[left-1][2],vtn[left-1][3],vtn[back-1][3],np+1,
						                                    vtn[left-1][6],vtn[left-1][7],vtn[back-1][7],np+2}));
						np+=2;
						
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 7 :
					 {
                        std::cout << i << " " << j << " " << k << std::endl;
                        D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(D[left-1][4]),-(D[back-1][3]),
						                                      nf+1,nf+2,nf+3}));
						nf+=3;
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[left-1][6],vte[bottom-1][10],vte[left-1][7],vte[bottom-1][11],
						                                    vte[back-1][7],ne+1,vte[left-1][11],vte[back-1][10],ne+2,ne+3}));
						ne+=3;
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[left-1][6],vtn[left-1][7],vtn[back-1][7],np+1}));
						np+=1;
						
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
				  }

				  tripletList.push_back(U(i,j,nv));
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
	  
	  //Now we build the rest of incidence matrices
	  Dt.resize(nf);
	  C.resize(nf);
	  Ct.resize(ne);
	  G.resize(ne);
	  Gt.resize(np);
	  M_h.resize(ne);
	  M_ni.resize(nf);
	  
	  for (int32_t i=0; i<nv; i++)
      {
         Dt[abs(D[i][0])].push_back(-i); 
		 Dt[abs(D[i][1])].push_back(-i);
         Dt[abs(D[i][2])].push_back(-i);
		 Dt[abs(D[i][3])].push_back( i);
		 Dt[abs(D[i][4])].push_back( i);
		 Dt[abs(D[i][5])].push_back( i);
		 
		 if (!C[abs(D[i][0])].size())
		 {
            C[abs(D[i][0])] = std::vector<int32_t>({vte[i][0],-vte[i][1],vte[i][3],-vte[i][5]});
			Ct[vte[i][0]].push_back( abs(D[i][0]));
			Ct[vte[i][1]].push_back(-abs(D[i][0]));
			Ct[vte[i][3]].push_back( abs(D[i][0]));
			Ct[vte[i][5]].push_back(-abs(D[i][0]));
		 }
		 if (!C[abs(D[i][1])].size())
		 {
            C[abs(D[i][1])] = std::vector<int32_t>({-vte[i][0],vte[i][2],-vte[i][4],vte[i][8]});
			Ct[vte[i][0]].push_back(-abs(D[i][1]));
			Ct[vte[i][2]].push_back( abs(D[i][1]));
			Ct[vte[i][4]].push_back(-abs(D[i][1]));
			Ct[vte[i][8]].push_back( abs(D[i][1]));
		 }
		 if (!C[abs(D[i][2])].size())
		 {
            C[abs(D[i][2])] = std::vector<int32_t>({vte[i][1],-vte[i][2],vte[i][6],-vte[i][9]});
			Ct[vte[i][1]].push_back( abs(D[i][2]));
			Ct[vte[i][2]].push_back(-abs(D[i][2]));
			Ct[vte[i][6]].push_back( abs(D[i][2]));
			Ct[vte[i][9]].push_back(-abs(D[i][2]));
		 }
		 if (!C[abs(D[i][3])].size())
		 {
            C[abs(D[i][3])] = std::vector<int32_t>({vte[i][3],-vte[i][4],vte[i][7],-vte[i][10]});
			Ct[vte[i][3]].push_back( abs(D[i][3]));
			Ct[vte[i][4]].push_back(-abs(D[i][3]));
			Ct[vte[i][7]].push_back( abs(D[i][3]));
			Ct[vte[i][10]].push_back(-abs(D[i][3]));
		 }
		 if (!C[abs(D[i][4])].size())
		 {
            C[abs(D[i][4])] = std::vector<int32_t>({-vte[i][5],vte[i][6],-vte[i][7],vte[i][11]});
			Ct[vte[i][5]].push_back(-abs(D[i][4]));
			Ct[vte[i][6]].push_back( abs(D[i][4]));
			Ct[vte[i][7]].push_back(-abs(D[i][4]));
			Ct[vte[i][11]].push_back(abs(D[i][4]));
		 }
		 if (!C[abs(D[i][5])].size())
		 {
            C[abs(D[i][5])] = std::vector<int32_t>({vte[i][8],-vte[i][9],vte[i][10],-vte[i][11]});
			Ct[vte[i][8]].push_back( abs(D[i][5]));
			Ct[vte[i][9]].push_back(-abs(D[i][5]));
			Ct[vte[i][10]].push_back( abs(D[i][5]));
			Ct[vte[i][11]].push_back(-abs(D[i][5]));
		 }
		 
		 if (!G[vte[i][0]].size())
		 {
            G[vte[i][0]]  = std::vector<int32_t>({-vtn[i][0],vtn[i][1]});
			Gt[vtn[i][0]].push_back(-vte[i][0]);
			Gt[vtn[i][1]].push_back( vte[i][0]);
		 }
		 if (!G[vte[i][1]].size())
		 {
            G[vte[i][1]]  = std::vector<int32_t>({-vtn[i][0],vtn[i][2]});
			Gt[vtn[i][0]].push_back(-vte[i][1]);
			Gt[vtn[i][2]].push_back( vte[i][1]);
		 }
		 if (!G[vte[i][2]].size())
		 {
            G[vte[i][0]]  = std::vector<int32_t>({-vtn[i][0],vtn[i][4]});
			Gt[vtn[i][0]].push_back(-vte[i][0]);
			Gt[vtn[i][4]].push_back( vte[i][0]);
		 }
		 if (!G[vte[i][3]].size())
		 {
            G[vte[i][3]]  = std::vector<int32_t>({-vtn[i][1],vtn[i][3]});
			Gt[vtn[i][1]].push_back(-vte[i][3]);
			Gt[vtn[i][3]].push_back( vte[i][3]);
		 }
		 if (!G[vte[i][4]].size())
		 {
            G[vte[i][4]]  = std::vector<int32_t>({-vtn[i][1],vtn[i][5]});
			Gt[vtn[i][1]].push_back(-vte[i][4]);
			Gt[vtn[i][5]].push_back( vte[i][4]);
		 }
		 if (!G[vte[i][5]].size())
		 {
            G[vte[i][5]]  = std::vector<int32_t>({-vtn[i][2],vtn[i][3]});
			Gt[vtn[i][2]].push_back(-vte[i][5]);
			Gt[vtn[i][3]].push_back( vte[i][5]);
		 }
		 if (!G[vte[i][6]].size())
		 {
            G[vte[i][6]]  = std::vector<int32_t>({-vtn[i][2],vtn[i][6]});
			Gt[vtn[i][2]].push_back(-vte[i][6]);
			Gt[vtn[i][6]].push_back( vte[i][6]);
		 }
		 if (!G[vte[i][7]].size())
		 {
            G[vte[i][7]]  = std::vector<int32_t>({-vtn[i][3],vtn[i][7]});
			Gt[vtn[i][3]].push_back(-vte[i][7]);
			Gt[vtn[i][7]].push_back( vte[i][7]);
		 }
		 if (!G[vte[i][8]].size())
		 {
            G[vte[i][8]]  = std::vector<int32_t>({-vtn[i][4],vtn[i][5]});
			Gt[vtn[i][4]].push_back(-vte[i][8]);
			Gt[vtn[i][5]].push_back( vte[i][8]);
		 }
		 if (!G[vte[i][9]].size())
		 {
            G[vte[i][9]]  = std::vector<int32_t>({-vtn[i][4],vtn[i][6]});
			Gt[vtn[i][4]].push_back(-vte[i][9]);
			Gt[vtn[i][6]].push_back( vte[i][9]);
		 }
		 if (!G[vte[i][10]].size())
		 {
            G[vte[i][10]]  = std::vector<int32_t>({-vtn[i][5],vtn[i][7]});
			Gt[vtn[i][5]].push_back(-vte[i][10]);
			Gt[vtn[i][7]].push_back( vte[i][10]);
		 }
		 if (!G[vte[i][11]].size())
		 {
            G[vte[i][11]]  = std::vector<int32_t>({-vtn[i][6],vtn[i][7]});
			Gt[vtn[i][6]].push_back(-vte[i][11]);
			Gt[vtn[i][7]].push_back( vte[i][11]);
		 }
      }
	  
	  std::cout << "ci arrivo qui?" << std::endl;
	  for (uint32_t i=0; i<nf; i++)
         M_ni[i]=avg_mu(i)/face_area(i);
	  
	  for (uint32_t i=0; i<ne; i++)
         M_h[i]=edge_len(i)/avg_eps(i);
   }
   
   T avg_mu(uint32_t index)
   {
      
   }

   T avg_eps(uint32_t index)
   {
      
   }
   
   bool Run(double simulation_time)
   {
      for (size_t i=0; i*t_step <= simulation_time; i++)
	  {
         for (size_t j=0; j<u.size(); j++)
		 {
			T i=0;
			for(size_t k=0; k < Ct[j].size(); k++)
				i+= Ct[j][k]/abs(Ct[j][k])*f[abs(Ct[j][k])-1];
            u[j]=u[j] + t_step*M_h[j]*i;
		 }
		 
         for (size_t j=0; j<f.size(); j++)
		 {
			T b = C[j][0]/abs(C[j][0])*u[abs(C[j][0])-1]+C[j][1]/abs(C[j][1])*u[abs(C[j][1])-1]+
			      C[j][2]/abs(C[j][2])*u[abs(C[j][2])-1]+C[j][3]/abs(C[j][3])*u[abs(C[j][3])-1];
            f[j]=f[j] - t_step*M_ni[j]*b;
		 }
	  }
   }

   uint32_t Volumes_size() { return D.size(); }
   uint32_t Surfaces_size() { return C.size(); }
   uint32_t Edges_size() { return G.size(); }
   uint32_t Points_size() { return pts.size(); }
   
   private:
   T xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz,L;
   std::vector<T> M_ni, M_h, u, f;
   double t_step;
   std::vector<std::vector<int32_t>> D,Dt,C,Ct,G,Gt;
   std::vector<Eigen::Vector3d> pts;
   
   double edge_len(uint32_t index)
   {
      return (pts[abs(G[index][0])]-pts[abs(G[index][1])]).norm();
   }
   
   double face_area(uint32_t index)
   {
      return edge_len(C[index][0])*edge_len(C[index][1]);
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
      if (!( pow(x,2)+pow(y,2)+pow(z,2) < 0.25 ))
         return false;
      return true;
   }
};
#endif
