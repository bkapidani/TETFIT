      //other instructions in constructor
	  // std::cout << "Qui?" << std::endl;
	  auto px = xmin;
	  auto py = ymin;
	  auto pz = zmin;  
	  
	  Lx = Ly = Lz = L;
	  t_step = 0.5*sqrt(pow(Lx,2)+pow(Ly,2)+pow(Lz,2))/c0/sqrt(3);
	  volume = Lx*Ly*Lz;
	  
	  uint32_t Nx = (fabs(xmax-xmin)) / Lx;// + 1;
	  uint32_t Ny = (fabs(ymax-ymin)) / Ly;// + 1;
	  uint32_t Nz = (fabs(zmax-zmin)) / Lz;// + 1;
	  
	  typedef Eigen::Triplet<uint32_t> Q;
      
	  Eigen::SparseMatrix<uint32_t> previous_layer(Nx,Ny);
	  
	  // previous_layer.SetZero();
	
      Eigen::Vector3d inc_x(Lx,0,0), inc_y(0,Ly,0), inc_z(0,0,Lz), dummy_vec;
	  bool not_found = true;
	  
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
							face_bars.push_back(dummy_vec);
						}
						
						vte.push_back(std::vector<int32_t>({ne+1,ne+2,ne+3,ne+4,ne+5,ne+6,
						                                    ne+7,ne+8,ne+9,ne+10,ne+11,ne+12}));
						ne+=12;
						
						for (uint32_t cnt=0; cnt<12; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						
						vtn.push_back(std::vector<int32_t>({np+1,np+2,np+3,np+4,
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],ne+1,vte[bottom-1][10],
						                                    ne+2,vte[bottom-1][11],ne+3,ne+4,ne+5,ne+6,ne+7,ne+8}));
						ne+=8;
						
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}

						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({vte[left-1][5],ne+1,vte[left-1][6],ne+2,vte[left-1][7],ne+3,
						                                    ne+4,ne+5,vte[left-1][11],ne+6,ne+7,ne+8}));
						ne+=8;
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						vtn.push_back(std::vector<int32_t>({vtn[left-1][2],vtn[left-1][3],np+1,np+2,
						                                    vtn[left-1][6],vtn[left-1][7],np+3,np+4}));
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[left-1][6],vte[bottom-1][10],vte[left-1][7],vte[bottom-1][11],
						                                    ne+1,ne+2,vte[left-1][11],ne+3,ne+4,ne+5}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[left-1][6],vtn[left-1][7],np+1,np+2}));
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({ne+1,vte[back-1][3],vte[back-1][4],ne+2,ne+3,ne+4,
						                                    vte[back-1][7],ne+5,ne+6,vte[back-1][10],ne+7,ne+8}));
						ne+=8;
						for (uint32_t cnt=0; cnt<8; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						vtn.push_back(std::vector<int32_t>({vtn[back-1][1],np+1,vtn[back-1][3],np+2,
						                                    vtn[back-1][5],np+3,vtn[back-1][7],np+4}));
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[back-1][4],vte[bottom-1][10],ne+1,vte[bottom-1][11],
						                                    vte[back-1][7],ne+2,ne+3,vte[back-1][10],ne+4,ne+5}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[back-1][5],np+1,vtn[back-1][7],np+2}));
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({vte[left-1][5],vte[back-1][3],vte[left-1][6],ne+1,vte[left-1][7],ne+2,
						                                    vte[back-1][7],ne+3,vte[left-1][11],vte[back-1][10],ne+4,ne+5}));
						ne+=5;
						for (uint32_t cnt=0; cnt<5; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						vtn.push_back(std::vector<int32_t>({vtn[left-1][2],vtn[left-1][3],vtn[back-1][3],np+1,
						                                    vtn[left-1][6],vtn[left-1][7],vtn[back-1][7],np+2}));
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
							face_bars.push_back(dummy_vec);
						}
						vte.push_back(std::vector<int32_t>({vte[bottom-1][8],vte[bottom-1][9],vte[left-1][6],vte[bottom-1][10],vte[left-1][7],vte[bottom-1][11],
						                                    vte[back-1][7],ne+1,vte[left-1][11],vte[back-1][10],ne+2,ne+3}));
						ne+=3;
						for (uint32_t cnt=0; cnt<3; cnt++)
						{
							Ct.push_back(dummyf);
							G.push_back(dummyf);
						}
						vtn.push_back(std::vector<int32_t>({vtn[bottom-1][4],vtn[bottom-1][5],vtn[bottom-1][6],vtn[bottom-1][7],
						                                    vtn[left-1][6],vtn[left-1][7],vtn[back-1][7],np+1}));
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
				  std::vector<Eigen::Vector3d> to_be_averaged;
				  for (auto nn : vtn[nv-1])
                     to_be_averaged.push_back(pts[abs(nn)-1]);
				  dual_pts.push_back(VectorAverage(to_be_averaged));
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
	  
	  std::vector<uint32_t> p_queue;
	  p_queue.push_back(1);
	  
	  //Now we build the rest of incidence matrices
	  // Dt.resize(nf);
	  // C.resize(nf);
	  // Ct.resize(ne);
	  // G.resize(ne);
	  // Gt.resize(np);
	  // M_h.resize(ne);
	  // M_ni.resize(nf);
	  // face_bars.resize(nf);
	  

	  for (int32_t i=0; i<nv; i++)
      {
	     // std::thread G_thread([&] {
		 if (!G[vte[i][0]-1].size())
		 {
            G[vte[i][0]-1]  = std::vector<int32_t>({-vtn[i][0],vtn[i][1]});
			Gt[vtn[i][0]-1].push_back(-vte[i][0]);
			Gt[vtn[i][1]-1].push_back( vte[i][0]);
		 }
		 if (!G[vte[i][1]-1].size())
		 {
            G[vte[i][1]-1]  = std::vector<int32_t>({-vtn[i][0],vtn[i][2]});
			Gt[vtn[i][0]-1].push_back(-vte[i][1]);
			Gt[vtn[i][2]-1].push_back( vte[i][1]);
		 }
		 if (!G[vte[i][2]-1].size())
		 {
            G[vte[i][2]-1]  = std::vector<int32_t>({-vtn[i][0],vtn[i][4]});
			Gt[vtn[i][0]-1].push_back(-vte[i][2]);
			Gt[vtn[i][4]-1].push_back( vte[i][2]);
		 }
		 if (!G[vte[i][3]-1].size())
		 {
            G[vte[i][3]-1]  = std::vector<int32_t>({-vtn[i][1],vtn[i][3]});
			Gt[vtn[i][1]-1].push_back(-vte[i][3]);
			Gt[vtn[i][3]-1].push_back( vte[i][3]);
		 }
		 if (!G[vte[i][4]-1].size())
		 {
            G[vte[i][4]-1]  = std::vector<int32_t>({-vtn[i][1],vtn[i][5]});
			Gt[vtn[i][1]-1].push_back(-vte[i][4]);
			Gt[vtn[i][5]-1].push_back( vte[i][4]);
		 }
		 if (!G[vte[i][5]-1].size())
		 {
            G[vte[i][5]-1]  = std::vector<int32_t>({-vtn[i][2],vtn[i][3]});
			Gt[vtn[i][2]-1].push_back(-vte[i][5]);
			Gt[vtn[i][3]-1].push_back( vte[i][5]);
		 }
		 if (!G[vte[i][6]-1].size())
		 {
            G[vte[i][6]-1]  = std::vector<int32_t>({-vtn[i][2],vtn[i][6]});
			Gt[vtn[i][2]-1].push_back(-vte[i][6]);
			Gt[vtn[i][6]-1].push_back( vte[i][6]);
		 }
		 if (!G[vte[i][7]-1].size())
		 {
            G[vte[i][7]-1]  = std::vector<int32_t>({-vtn[i][3],vtn[i][7]});
			Gt[vtn[i][3]-1].push_back(-vte[i][7]);
			Gt[vtn[i][7]-1].push_back( vte[i][7]);
		 }
		 if (!G[vte[i][8]-1].size())
		 {
            G[vte[i][8]-1]  = std::vector<int32_t>({-vtn[i][4],vtn[i][5]});
			Gt[vtn[i][4]-1].push_back(-vte[i][8]);
			Gt[vtn[i][5]-1].push_back( vte[i][8]);
		 }
		 if (!G[vte[i][9]-1].size())
		 {
            G[vte[i][9]-1]  = std::vector<int32_t>({-vtn[i][4],vtn[i][6]});
			Gt[vtn[i][4]-1].push_back(-vte[i][9]);
			Gt[vtn[i][6]-1].push_back( vte[i][9]);
		 }
		 if (!G[vte[i][10]-1].size())
		 {
            G[vte[i][10]-1]  = std::vector<int32_t>({-vtn[i][5],vtn[i][7]});
			Gt[vtn[i][5]-1].push_back(-vte[i][10]);
			Gt[vtn[i][7]-1].push_back( vte[i][10]);
		 }
		 if (!G[vte[i][11]-1].size())
		 {
            G[vte[i][11]-1]  = std::vector<int32_t>({-vtn[i][6],vtn[i][7]});
			Gt[vtn[i][6]-1].push_back(-vte[i][11]);
			Gt[vtn[i][7]-1].push_back( vte[i][11]);
		 }
         // });		 
         // std::thread C_thread([&] {
		 if (!C[abs(D[i][0])-1].size())
		 {
            C[abs(D[i][0])-1] = std::vector<int32_t>({vte[i][0],-vte[i][1],vte[i][3],-vte[i][5]});
			face_bars[abs(D[i][0])-1] = face_barycenter(abs(D[i][0])-1);
			Ct[vte[i][0]-1].push_back( abs(D[i][0]));
			Ct[vte[i][1]-1].push_back(-abs(D[i][0]));
			Ct[vte[i][3]-1].push_back( abs(D[i][0]));
			Ct[vte[i][5]-1].push_back(-abs(D[i][0]));
		 }
		 if (!C[abs(D[i][1])-1].size())
		 {
            C[abs(D[i][1])-1] = std::vector<int32_t>({-vte[i][0],vte[i][2],-vte[i][4],vte[i][8]});
			face_bars[abs(D[i][1])-1] = face_barycenter(abs(D[i][1])-1);
			Ct[vte[i][0]-1].push_back(-abs(D[i][1]));
			Ct[vte[i][2]-1].push_back( abs(D[i][1]));
			Ct[vte[i][4]-1].push_back(-abs(D[i][1]));
			Ct[vte[i][8]-1].push_back( abs(D[i][1]));
		 }
		 if (!C[abs(D[i][2])-1].size())
		 {
            C[abs(D[i][2])-1] = std::vector<int32_t>({vte[i][1],-vte[i][2],vte[i][6],-vte[i][9]});
			face_bars[abs(D[i][2])-1] = face_barycenter(abs(D[i][2])-1);
			Ct[vte[i][1]-1].push_back( abs(D[i][2]));
			Ct[vte[i][2]-1].push_back(-abs(D[i][2]));
			Ct[vte[i][6]-1].push_back( abs(D[i][2]));
			Ct[vte[i][9]-1].push_back(-abs(D[i][2]));
		 }
		 if (!C[abs(D[i][3])-1].size())
		 {
            C[abs(D[i][3])-1] = std::vector<int32_t>({vte[i][3],-vte[i][4],vte[i][7],-vte[i][10]});
			face_bars[abs(D[i][3])-1] = face_barycenter(abs(D[i][3])-1);
			Ct[vte[i][3]-1].push_back( abs(D[i][3]));
			Ct[vte[i][4]-1].push_back(-abs(D[i][3]));
			Ct[vte[i][7]-1].push_back( abs(D[i][3]));
			Ct[vte[i][10]-1].push_back(-abs(D[i][3]));
		 }
		 if (!C[abs(D[i][4])-1].size())
		 {
            C[abs(D[i][4])-1] = std::vector<int32_t>({-vte[i][5],vte[i][6],-vte[i][7],vte[i][11]});
			face_bars[abs(D[i][4])-1] = face_barycenter(abs(D[i][4])-1);
			Ct[vte[i][5]-1].push_back(-abs(D[i][4]));
			Ct[vte[i][6]-1].push_back( abs(D[i][4]));
			Ct[vte[i][7]-1].push_back(-abs(D[i][4]));
			Ct[vte[i][11]-1].push_back(abs(D[i][4]));
		 }
		 if (!C[abs(D[i][5])-1].size())
		 {
            C[abs(D[i][5])-1] = std::vector<int32_t>({vte[i][8],-vte[i][9],vte[i][10],-vte[i][11]});
			face_bars[abs(D[i][5])-1] = face_barycenter(abs(D[i][5])-1);
			Ct[vte[i][8]-1].push_back( abs(D[i][5]));
			Ct[vte[i][9]-1].push_back(-abs(D[i][5]));
			Ct[vte[i][10]-1].push_back( abs(D[i][5]));
			Ct[vte[i][11]-1].push_back(-abs(D[i][5]));
		 }
		 // });
	  

	     // G_thread.join();
	     // C_thread.join();
		 
		 
      }
	  
	  // for (size_t k=0; k<nf; k++)
	  // {
         // assert(Dt[k].size()>=1 && Dt[k].size()<=2);
		 // assert(C[k].size() == 4);
      // }
	  // for (size_t k=0; k<ne; k++)
	  // {
         // assert(Ct[k].size()>=2 && Ct[k].size()<=4);
		 // assert(G[k].size() == 2);
      // }
	  // for (size_t k=0; k<np; k++)
	  // {
         // assert(Gt[k].size()>=3 && Gt[k].size()<=6);
      // }
	  
      t_mesh.toc();
	  std::cout << "elements: " << nv << " surfaces: " << nf << " edges: " << ne << " nodes: " << np << " dofs: " << nf+ne << std::endl;
	  std::cout << "Meshing takes " << t_mesh << " seconds!" << std::endl;
	  // std::cout << pts[pts.size()-1] << std::endl;
	  
	  t_mesh.tic();
	  for (uint32_t i=0; i<nf; i++)
	  {
         M_ni.push_back(avg_ni(i)/face_area(i));
		 F.push_back(0);
	  }
	  for (uint32_t i=0; i<ne; i++)
	  {
         M_h.push_back(edge_len(i)/avg_eps(i));
		 
		 U.push_back(0);
		 switch (Ct[i].size())
		 {
			case 2 :
			{
				bc.push_back(0);
				is_boundary.push_back(true);
				break;
			}
			case 3 :
			{
				is_boundary.push_back(true);
				auto bar = edge_barycenter(i);
				if (bar(2) <= zmin)
					bc.push_back(sin( bar(0)*pi/0.05 )*edge_vector(i)(1));
				else
					bc.push_back(0);

				break;
			}
			case 4 :
			{
				bc.push_back(1);
				is_boundary.push_back(false);
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
   
   T avg_ni(uint32_t index)
   {
	  T ret=0;
	  for (auto v : Dt[index])
	  {
         auto label = abs(v) - 1;
		 ret += (dual_pts[label]-face_bars[index]).norm()/mu[material[label]];
	  }
      return ret;
   }

   T avg_eps(uint32_t index)
   {
	  T ret=0;
	  std::map<uint32_t,std::vector<uint32_t>> dual_faces;
	  for (auto f : Ct[index])
	  {
		 auto abs_f = abs(f)-1;
		 for (auto v : Dt[abs_f])
            dual_faces[abs(v)-1].push_back(abs_f);
	  }
	  for (auto v : dual_faces)
	  {
         auto label = v.first;
		 ret += (dual_pts[label]-face_bars[v.second[0]]).norm()*
		        (dual_pts[label]-face_bars[v.second[1]]).norm()*
				 epsilon[material[label]];
	  }
	  
      return ret;
   }
   
   bool Run(double simulation_time)
   {
	  double step_time_average=0;
	  timecounter step_cost;
	  size_t i;
	  std::vector<T> numeric_values,numeric_times;
      for (i=0; i*t_step <= simulation_time; i++)
	  {
		 step_cost.tic();
		 T time_function=sin(2*pi*freq*i*t_step);
         for (size_t j=0; j<U.size(); j++)
		 {
			
			if (is_boundary[j])
			{
				U[j]=time_function*bc[j];
				// std::cout << "U[" << j << "] = " << U[j] << std::endl;
			}
			else
			{
               // T Dpartial=0;
			   // for(size_t k=0; k < 4; k++)
			      // Dpartial+=   Ct[j][k]/abs(Ct[j][k])*F[abs(Ct[j][k])-1];
			   T Dpartial = sgn(Ct[j][0])*F[abs(Ct[j][0])-1]+sgn(Ct[j][1])*F[abs(Ct[j][1])-1]+
			                sgn(Ct[j][2])*F[abs(Ct[j][2])-1]+sgn(Ct[j][3])*F[abs(Ct[j][3])-1];
               U[j]=U[j] + t_step*M_h[j]*Dpartial;
			   
			}
			// std::cout << "U[" << j << "] = " << U[j] << " ";
		 }
		 
		 // std::cout << std::endl;
         for (size_t j=0; j<F.size(); j++)
		 {
			// T Bpartial = 0;
			// for(size_t k=0; k < 4; k++)
			   // Bpartial += C[j][k]/abs(C[j][k])*U[abs(C[j][k])-1];
			T Bpartial   = sgn(C[j][0])*U[abs(C[j][0])-1]+sgn(C[j][1])*U[abs(C[j][1])-1]+
			               sgn(C[j][2])*U[abs(C[j][2])-1]+sgn(C[j][3])*U[abs(C[j][3])-1];
            F[j]= F[j] - t_step*M_ni[j]*Bpartial;
		 }
		 auto num_val = GetElectricfield(probe_elem);
		 numeric_values.push_back(num_val(1));
		 numeric_times.push_back(i*t_step);
		 step_cost.toc();
		 step_time_average += (duration_cast<duration<double>>(step_cost.elapsed())).count();
	  }
	  std::ofstream os;
	  os.open("numeric_FIT.dat");
	  for (size_t k=0; k < numeric_values.size(); k++)
		  os << numeric_times[k] << " " << numeric_values[k] << std::endl;
	  os.close();
	  std::cout << "Time step takes (average) " << step_time_average/(double(i)) << " seconds" << std::endl;
	  std::cout << "Total running time is "     << step_time_average << " seconds" << std::endl;
   }
   
   uint32_t Volumes_size() { return D.size(); }
   uint32_t Surfaces_size() { return C.size(); }
   uint32_t Edges_size() { return G.size(); }
   uint32_t Points_size() { return pts.size(); }
   
   private:
   T xmin,xmax,ymin,ymax,zmin,zmax,Lx,Ly,Lz,L, volume;
   uint32_t probe_elem;
   std::vector<T> bc;
   std::vector<bool> is_boundary;
   std::vector<T> M_ni, M_h, U, F;
   std::map<uint32_t,double> epsilon,mu;
   double t_step, freq;
   std::vector<uint32_t> material;
   std::vector<std::vector<int32_t>> D,Dt,C,Ct,G,Gt;
   std::vector<Eigen::Vector3d> pts, dual_pts, face_bars;
   std::vector<std::vector<int32_t>> vte,vtn;
   
   int32_t sgn(int32_t val)
   {
      return (0 < val) - (val < 0);
   }
   
   Eigen::Vector3d GetElectricfield(uint32_t cube)
   {
	  // auto volume = edge_len(vtn[cube][0]-1)*edge_len(vtn[cube][1]-1)*edge_len(vtn[cube][2]-1);
	  std::vector<Eigen::Vector3d> edge_mp,face_mp;
	  std::vector<T> u;
	  
	  for (uint32_t i=0; i<12; i++)
	  {
		  edge_mp.push_back(edge_barycenter(vte[cube][i]-1));
		  u.push_back(U[vte[cube][i]-1]);
	  }
	  for (uint32_t i=0; i<6; i++)
		  face_mp.push_back(face_bars[abs(D[cube][i])-1]);
	  
      Eigen::Vector3d ret = (u[0]*(face_mp[0]-edge_mp[0]).cross(face_mp[1]-edge_mp[0])    +
	                         u[1]*(face_mp[2]-edge_mp[1]).cross(face_mp[0]-edge_mp[1])    +
	                         u[2]*(face_mp[1]-edge_mp[2]).cross(face_mp[2]-edge_mp[2])    +
	                         u[3]*(face_mp[0]-edge_mp[3]).cross(face_mp[3]-edge_mp[3])    +
	                         u[4]*(face_mp[3]-edge_mp[4]).cross(face_mp[1]-edge_mp[4])    +
	                         u[5]*(face_mp[4]-edge_mp[5]).cross(face_mp[0]-edge_mp[5])    +
	                         u[6]*(face_mp[2]-edge_mp[6]).cross(face_mp[4]-edge_mp[6])    +
	                         u[7]*(face_mp[4]-edge_mp[7]).cross(face_mp[3]-edge_mp[7])    +
	                         u[8]*(face_mp[1]-edge_mp[8]).cross(face_mp[5]-edge_mp[8])    +
	                         u[9]*(face_mp[5]-edge_mp[9]).cross(face_mp[2]-edge_mp[9])    +
	                         u[10]*(face_mp[3]-edge_mp[10]).cross(face_mp[5]-edge_mp[10]) +
	                         u[11]*(face_mp[5]-edge_mp[11]).cross(face_mp[4]-edge_mp[11]))/volume;
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
   
   double edge_len(uint32_t index)
   {
	  // std::cout << "index: " << index << std::endl;
	  // std::cout << abs(G[index][0])-1 << " " << abs(G[index][1])-1 << std::endl;
      return (pts[abs(G[index][0])-1]-pts[abs(G[index][1])-1]).norm();
	  
   }
   
   double face_area(uint32_t index)
   {
	  // std::cout << abs(C[index][0])-1 << " " << abs(C[index][1])-1 << " " << abs(C[index][2])-1 << " " << abs(C[index][3])-1 << std::endl;
      return edge_len(abs(C[index][0])-1)*edge_len(abs(C[index][1])-1);
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
