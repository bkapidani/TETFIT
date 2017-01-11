      //other instructions in constructor
	  auto px = xmin;
	  auto py = ymin;
	  auto pz = zmin;  
	  
	  uint32_t Nx = (fabs(xmax-xmin)) / L;
	  uint32_t Ny = (fabs(ymax-ymin)) / L;
	  uint32_t Nz = (fabs(zmax-zmin)) / L;
	  
	  typedef Eigen::Triplet<uint32_t> U;
      std::vector<U> tripletList;
      tripletList.reserve(Nx*Ny);
      
	  Eigen::SparseMatrix<uint32_t> previous_layer(Nx,Ny);
	  
	  // previous_layer.SetZero();
	
      Eigen::Vector3d inc_x(L,0,0);
	  Eigen::Vector3d inc_y(0,L,0);
	  Eigen::Vector3d inc_z(0,0,L);
	  std::vector<std::vector<uint32_t>> vte,vtn;
	  
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
                        D.push_back(std::vector<int32_t>({-(nf+1),-(nf+2),-(nf+3),
						                                      nf+4,nf+5,nf+6}));
						nf+=6;
						vte.push_back(std::vector<uint32_t>({ne+1,ne+2,ne+3,ne+4,ne+5,ne+6,
						                                    ne+7,ne+8,ne+9,ne+10,ne+11,ne+12}));
						ne+=12;
						vtn.push_back(std::vector<uint32_t>({np+1,np+2,np+3,np+4,
						                                    np+5,np+6,np+7,np+8}));
						np+=8;
						
                        pts.push_back(pp);						//nte.push_back(dummyf);
						pts.push_back(pp+inc_x);				//nte.push_back(dummyf);
						pts.push_back(pp+inc_y);				//nte.push_back(dummyf);
						pts.push_back(pp+inc_x+inc_y);			//nte.push_back(dummyf);
                        pts.push_back(pp+inc_z);				//nte.push_back(dummyf);
						pts.push_back(pp+inc_z+inc_x);			//nte.push_back(dummyf);
						pts.push_back(pp+inc_z+inc_y);			//nte.push_back(dummyf);
						pts.push_back(pp+inc_z+inc_x+inc_y);	//nte.push_back(dummyf);
						
						break;
					 }
					 case 1 :
					 {
						D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(nf+1),-(nf+2),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						vte.push_back(std::vector<uint32_t>({vte[bottom-1][8],vte[bottom-1][9],ne+1,vte[bottom-1][10],
						                                    ne+2,vte[bottom-1][11],ne+3,ne+4,ne+5,ne+6,ne+7,ne+8}));
						ne+=8;
						vtn.push_back(std::vector<uint32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
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
						D.push_back(std::vector<int32_t>({-(nf+1),-(D[left-1][5]),-(nf+2),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						vte.push_back(std::vector<uint32_t>({vte[left-1][5],ne+1,vte[left-1][6],ne+2,vte[left-1][7],ne+3,
						                                    ne+4,ne+5,vte[left-1][11],ne+6,ne+7,ne+8}));
						ne+=8;
						vtn.push_back(std::vector<uint32_t>({vtn[left-1][2],vtn[left-1][3],np+1,np+2,
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
                        D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(D[left-1][4]),-(nf+1),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						vte.push_back(std::vector<uint32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[left-1][6],vte[bottom-1][10],vte[left-1][7],vte[bottom-1][11],
						                                    ne+1,ne+2,vte[left-1][11],ne+3,ne+4,ne+5}));
						ne+=5;
						vtn.push_back(std::vector<uint32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[left-1][6],vtn[left-1][7],np+1,np+2}));
						np+=2;
						 
						pts.push_back(pp+inc_z+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 4 :
					 {
						D.push_back(std::vector<int32_t>({-(nf+1),-(nf+2),-(D[back-1][3]),
						                                      nf+3,nf+4,nf+5}));
						nf+=5;
						vte.push_back(std::vector<uint32_t>({ne+1,vte[back-1][3],vte[back-1][4],ne+2,ne+3,ne+4,
						                                    vte[back-1][7],ne+5,ne+6,vte[back-1][10],ne+7,ne+8}));
						ne+=8;
						vtn.push_back(std::vector<uint32_t>({vtn[back-1][1],np+1,vtn[back-1][3],np+2,
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
						D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(nf+1),-(D[back-1][3]),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						vte.push_back(std::vector<uint32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[back-1][4],vte[bottom-1][10],ne+1,vte[bottom-1][11],
						                                    vte[back-1][7],ne+2,ne+3,vte[back-1][10],ne+4,ne+5}));
						ne+=5;
						vtn.push_back(std::vector<uint32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[back-1][5],np+1,vtn[back-1][7],np+2}));
						np+=2;
						
						pts.push_back(pp+inc_z+inc_x);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 6 :
					 {
						D.push_back(std::vector<int32_t>({-(nf+1),-(D[left-1][5]),-(D[back-1][3]),
						                                      nf+2,nf+3,nf+4}));
						nf+=4;
						vte.push_back(std::vector<uint32_t>({vte[left-1][5],vte[back-1][3],vte[left-1][6],ne+1,vte[left-1][7],ne+2,
						                                    vte[back-1][7],ne+3,vte[left-1][11],vte[back-1][10],ne+4,ne+5}));
						ne+=5;
						vtn.push_back(std::vector<uint32_t>({vtn[left-1][2],vtn[left-1][3],vtn[back-1][3],np+1,
						                                    vtn[left-1][6],vtn[left-1][7],vtn[back-1][7],np+2}));
						np+=2;
						
						pts.push_back(pp+inc_x+inc_y);
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
					 case 7 :
					 {
                        D.push_back(std::vector<int32_t>({-(D[bottom-1][5]),-(D[left-1][4]),-(D[back-1][3]),
						                                      nf+1,nf+2,nf+3}));
						nf+=3;
						vte.push_back(std::vector<uint32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[left-1][6],vte[bottom-1][10],vte[left-1][7],vte[bottom-1][11],
						                                    vte[back-1][7],ne+1,vte[left-1][11],vte[back-1][10],ne+2,ne+3}));
						ne+=3;
						vtn.push_back(std::vector<uint32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[left-1][6],vtn[left-1][7],vtn[back-1][7],np+1}));
						np+=1;
						
						pts.push_back(pp+inc_z+inc_x+inc_y);
						
						break;
					 }
				  }

				  tripletList.push_back(U(i,j,nv));
				  this_col[i]=nv;
			   }
			   
               px+=L;
			}
			
			old_col=std::move(this_col);
			py+=L;
		 }

         this_layer.setFromTriplets(tripletList.begin(), tripletList.end());
		 previous_layer=std::move(this_layer);
		 pz+=L;
      }
	  
	  //Now we build the rest of incidence matrices
	  Dt.resize(nf);
	  C.resize(nf);
	  Ct.resize(ne);
	  G.resize(ne);
	  Gt.resize(np);
	  
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
	  
	   
   }

   private:
   T xmin,xmax,ymin,ymax,zmin,zmax,L;
   std::vector<std::vector<int32_t>> D,Dt,C,Ct,G,Gt;
   std::vector<Eigen::Vector3d> pts;
   
   bool IsGridVol(T x, T y, T z)
   {
	   if (!IsGridPoint(x,y,z))
		   return false;
	   if (!IsGridPoint(x+L,y,z))
		   return false;
	   if (!IsGridPoint(x,y+L,z))
		   return false;
	   if (!IsGridPoint(x,y,z+L))
		   return false;
	   if (!IsGridPoint(x+L,y+L,z))
		   return false;
	   if (!IsGridPoint(x+L,y,z+L))
		   return false;
	   if (!IsGridPoint(x,y+L,z+L))
		   return false;
	   if (!IsGridPoint(x+L,y+L,z+L))
		   return false;
	   return true;
   }
   bool IsGridPoint(T x, T y, T z)
   {
      //check if it satisfies the set of inequalities
