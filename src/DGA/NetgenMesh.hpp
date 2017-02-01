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

//file NetgenMesh.hpp
#include <iostream>
#include <cstdlib>
#include <utility> 
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <map>
#include <set>
#include <limits>
#include <sstream>
#include <future>
#include <chrono>
#include <mutex>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "timecounter.h"
#include "sgnint32_t.hpp"
#include "mapped_file.h"
#include "strtot.hpp"

const long double PI = 3.141592653589793238L;

namespace parser {

template<typename T>
std::tuple<T, T, T>
read_point_line(const char *str, char **endptr)
{
    T t1, t2, t3;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    
    //return std::make_tuple(t1/1000.0, t2/1000.0, t3/1000.0);
    return std::make_tuple(t1, t2, t3);
}

template<typename T>
std::tuple<T, T, T, T, T>
read_tetrahedron_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4, t5;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    t5 = strtot<T>(*endptr, endptr);
    
    return std::make_tuple(t1, t2-1, t3-1, t4-1, t5-1);
}
    
template<typename T>
std::tuple<T, T, T, T>
read_triangle_line(const char *str, char **endptr)
{
    T t1, t2, t3, t4;
    
    t1 = strtot<T>(str, endptr);
    t2 = strtot<T>(*endptr, endptr);
    t3 = strtot<T>(*endptr, endptr);
    t4 = strtot<T>(*endptr, endptr);
    
    return std::make_tuple(t1, t2-1, t3-1, t4-1);
}


template<typename T>
void
sort_uniq(std::vector<T>& v)
{
    std::sort(v.begin(), v.end());
    auto uniq_iter = std::unique(v.begin(), v.end());
    v.erase(uniq_iter, v.end());
}
    
} //namespace priv

using cluster_list     		= std::vector<sgnint32_t<int32_t>>; 
using h1_2d_basis 	   		= std::vector<std::vector<int32_t>>;
using thinned_currents 		= std::vector<std::vector<int32_t>>;
using volume_type 			= std::tuple<uint32_t,uint32_t,uint32_t,uint32_t>;
using surface_type 			= std::tuple<uint32_t,uint32_t,uint32_t>;
using label_surface_type 	= std::pair<surface_type,uint32_t>;
using edge_type 			= std::tuple<uint32_t,uint32_t>;
using label_edge_type 		= std::pair<edge_type,uint32_t>;
using label_node_type 		= std::pair<uint32_t,uint32_t>;
using tm_tuple 				= std::tuple< volume_type, uint32_t >;
using sm_tuple 				= std::tuple< surface_type, uint32_t >;

class NetgenMesh
{	
	public:
	NetgenMesh(std::string meshfile)
	{
		timecounter t_read;
		t_read.tic();
		read_mesh(meshfile);
		t_read.toc();
		std::cout << "Loading complex took: " << t_read << " s" << std::endl;
	}

	void unique(std::vector<label_surface_type>& arr, std::vector<uint32_t>& new_labels)
	{
		if (!arr.size())
			throw std::invalid_argument("Array must be nonempty!");
		
		uint32_t left = 0;
		uint32_t right = arr.size()-1;

		struct {
			bool operator()(const label_surface_type& t1, const label_surface_type& t2)
			{
				return (t1.first < t2.first);
			}
		} surfcomp;
		
		std::sort(arr.begin(),arr.end(),surfcomp);	
		// std::vector<uint32_t> new_labels(arr.size(),0);

		// std::cout << std::endl;
		
		uint32_t itor = left;
		new_labels[arr[0].second]=itor;
		surfaces.push_back(arr[0].first);
		
		left++;
		while (left<=right)
		{
			if (arr[left].first == arr[left-1].first)
			{
				new_labels[arr[left].second]=itor;
			}
			else
			{
				itor++;
				new_labels[arr[left].second]=itor;
				surfaces.push_back(arr[left].first);
			}
			
			left++;
		}
		
		// arr = new_arr;
		// labels = new_labels;
		
		return;
		
	}

	void unique(std::vector<label_edge_type>& arr, std::vector<uint32_t>& new_labels)
	{
		if (!arr.size())
			throw std::invalid_argument("Array must be nonempty!");
		
		uint32_t left = 0;
		uint32_t right = arr.size()-1;

		struct {
			bool operator()(const label_edge_type& t1, const label_edge_type& t2)
			{
				return (t1.first < t2.first);
			}
		} edgecomp;
		
		std::sort(arr.begin(),arr.end(),edgecomp);	
		// std::vector<uint32_t> new_labels(arr.size(),0);
		
		uint32_t itor = left;
		new_labels[arr[0].second]=itor;
		edges.push_back(arr[0].first);
		
		left++;
		while (left<=right)
		{
			if (arr[left].first == arr[left-1].first)
			{
				new_labels[arr[left].second]=itor;
			}
			else
			{
				itor++;
				new_labels[arr[left].second]=itor;
				edges.push_back(arr[left].first);
			}
			
			left++;
		}
		
		// arr = new_arr;
		// labels = new_labels;
		
		return;
		
	}

	bool read_mesh(const std::string& _filename)
	{	
		timecounter tc, tctot;
		
		/* Open file */
		if (_filename.size() == 0)
		{
			std::cout << "Invalid mesh file name" << std::endl;
			return false;
		}
		
		uint32_t	lines, linecount;
		
		mapped_file mf(_filename);
		
		std::cout << " * * * Reading NETGEN format mesh * * * " << std::endl;
		
		tctot.tic();
		
		/************************ Read points ************************/
		linecount = 0;
		
		const char *data = mf.mem();
		char *endptr;
		
		lines = strtot<uint32_t>(data, &endptr);
		
		// pts.reserve(lines);
		
		tc.tic();
		while (linecount < lines)
		{
			if ( (linecount%100000) == 0 )
			{
				std::cout << "Reading points: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}

			auto t = parser::read_point_line<double>(endptr, &endptr);
			
			
			Eigen::Vector3d point = { std::get<0>(t),std::get<1>(t),std::get<2>(t) };		
			pts.push_back(point);
			
			/* Do something with that point */
			
			linecount++;
		}
		tc.toc();
		
		std::cout << "Reading points: " << linecount;
		std::cout << "/" << lines << " - " << tc << " seconds" << std::endl;
		
		/************************ Read tetrahedra ************************/
		linecount = 0;
		
		lines = strtot<uint32_t>(endptr, &endptr);
		std::vector< tm_tuple > temp_tet;
		temp_tet.reserve(lines);
		
		tc.tic();
		while (linecount < lines)
		{
			if ( (linecount%100000) == 0 )
			{
				std::cout << "Reading tetrahedra: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}
			
			auto t = parser::read_tetrahedron_line<uint32_t>(endptr, &endptr);
			
			//auto t = parser::read_array<uint32_t, 5>(mf.get_line());
			std::vector<uint32_t> vectet(4);
			vectet[0] = std::get<1>(t);
			vectet[1] = std::get<2>(t);
			vectet[2] = std::get<3>(t);
			vectet[3] = std::get<4>(t);
			std::sort(vectet.begin(),vectet.end());
			
			uint32_t       p0(vectet[0]);
			uint32_t       p1(vectet[1]);
			uint32_t       p2(vectet[2]);
			uint32_t       p3(vectet[3]);
			uint32_t       d(std::get<0>(t));
			
			auto tuple = std::make_tuple(volume_type(p0, p1, p2, p3), d);
			temp_tet.push_back( tuple );
			
			linecount++;
		}
		tc.toc();
		
		std::cout << "Reading tetrahedra: " << linecount;
		std::cout << "/" << lines  << " - " << tc << " seconds" << std::endl;
		
		/************************ Sort ************************/
		std::cout << "Sorting data...";
		std::cout.flush();
		
		tc.tic();
		
		/* sort tetrahedra, make unique and move them in geometry */
			
		struct {
			bool operator()(const tm_tuple& t1, const tm_tuple& t2)
			{
				return (std::get<0>(t1) < std::get<0>(t2));
			}
		} mycomp;

		std::sort(temp_tet.begin(), temp_tet.end(), mycomp);
		
		// _vtf_list.resize(lines);
		uint32_t tot=0;	

		std::vector<label_surface_type> temp_tri0;
		temp_tri0.resize(4*lines);
		std::vector<int32_t> vol_signs;
		
		// vol_signs.reserve(lines);
		// volumes.reserve(lines);
		// domains.reserve(lines);
		
		for (auto tet : temp_tet)
		{
			auto t = std::get<0>(tet);
			volumes.push_back(t);
			
			uint32_t       p0(std::get<0>(t));
			uint32_t       p1(std::get<1>(t));
			uint32_t       p2(std::get<2>(t));
			uint32_t       p3(std::get<3>(t));

			temp_tri0[tot]           = std::make_pair(surface_type(p0, p1, p2),tot);
			temp_tri0[tot+lines]     = std::make_pair(surface_type(p0, p1, p3),tot+lines);
			temp_tri0[tot+2*lines]   = std::make_pair(surface_type(p0, p2, p3),tot+2*lines);
			temp_tri0[tot+3*lines]   = std::make_pair(surface_type(p1, p2, p3),tot+3*lines);
			tot++;
		
			Eigen::Vector3d v1 = pts[p1] - pts[p0];
			Eigen::Vector3d v2 = pts[p2] - pts[p0];
			Eigen::Vector3d v3 = pts[p3] - pts[p0];		

			int32_t sgn  = v1.dot(v2.cross(v3))/double(6)>0? 1 : -1;
			vol_signs.push_back(sgn);
			domains.push_back(std::get<1>(tet));
		}
		std::vector<tm_tuple>().swap(temp_tet);
		std::vector<uint32_t> labels(4*lines);
		surfaces.reserve(4*lines);
		unique(temp_tri0, labels); //this also fills the surfaces vector
		std::vector<label_surface_type>().swap(temp_tri0);
		
		_ftv_list.resize(surfaces.size());
		// _vtf_list.reserve(volumes.size());
		
		for (size_t k=0; k<lines; k++)
		{
			// std::cout << labels[k] << " " << _ftv_list.size() << std::endl;
			sgnint32_t<int32_t> v1(k,-vol_signs[k]);
			sgnint32_t<int32_t> v2(k,vol_signs[k]);
			
			_ftv_list[labels[k]].push_back(v1); 
			_ftv_list[labels[k+lines]].push_back(v2); 
			_ftv_list[labels[k+2*lines]].push_back(v1);
			_ftv_list[labels[k+3*lines]].push_back(v2);		
			
			sgnint32_t<int32_t> f1(labels[k],-vol_signs[k]);
			sgnint32_t<int32_t> f2(labels[k+lines],vol_signs[k]);
			sgnint32_t<int32_t> f3(labels[k+2*lines],-vol_signs[k]);
			sgnint32_t<int32_t> f4(labels[k+3*lines],vol_signs[k]);
			
			std::vector<sgnint32_t<int32_t>> dummy(4);		
			_vtf_list.push_back(dummy);
			_vtf_list[k][0] = f1; 
			_vtf_list[k][1] = f2; 
			_vtf_list[k][2] = f3;
			_vtf_list[k][3] = f4;
		}
		
		std::vector<uint32_t>().swap(labels);
		lines = surfaces.size();
		// intersurface.resize(lines,0);
		// _fte_list.resize(lines);
		tot = 0;
		std::vector<label_edge_type> temp_edge0(3*lines);
		
		for (auto t : surfaces)
		{
			uint32_t       p0(std::get<0>(t));
			uint32_t       p1(std::get<1>(t));
			uint32_t       p2(std::get<2>(t));

			temp_edge0[tot]          =  std::make_pair(edge_type(p0, p1),tot);
			temp_edge0[tot+lines]    =  std::make_pair(edge_type(p0, p2),tot+lines);
			temp_edge0[tot+2*lines]  =  std::make_pair(edge_type(p1, p2),tot+2*lines);
			tot++;
		}
		
		std::vector<uint32_t> e_labels(3*lines);
		unique(temp_edge0, e_labels);
		std::vector<label_edge_type>().swap(temp_edge0);
		_etf_list.resize(edges.size());
		// _fte_list.reserve(surfaces.size());
		// _etn_list.resize(edges.size());
		// physical_edges.resize(edges_size(),0);
		
		for (size_t k=0; k<lines; k++)
		{
			
			sgnint32_t<int32_t> f1(k, 1);
			sgnint32_t<int32_t> f2(k,-1);
			sgnint32_t<int32_t> f3(k, 1);
			
			_etf_list[e_labels[k]].push_back(f1);
			_etf_list[e_labels[k+lines]].push_back(f2); 
			_etf_list[e_labels[k+2*lines]].push_back(f3);

			// std::cout << labels[k] << " " << labels[k+lines] << " " << labels[k+2*lines] << std::endl;
			
			sgnint32_t<int32_t> e1(e_labels[k],1);
			sgnint32_t<int32_t> e2(e_labels[k+lines],-1);
			sgnint32_t<int32_t> e3(e_labels[k+2*lines],1);
			
			std::vector<sgnint32_t<int32_t>> dummy(3);		
			_fte_list.push_back(dummy);
			
			_fte_list[k][0] = e1;
			_fte_list[k][1] = e2;
			_fte_list[k][2] = e3;
			
			auto vols = _ftv_list[k];
			
			switch (vols.size()) 
			{
				case 2: 
				{
					auto vol1= abs(*vols.begin());
					auto vol2= abs(*(std::prev(vols.end())));
					
					if (domains[vol1] != domains[vol2])
					{
						// intersurface[k]++;
						for (auto ee : _fte_list[k])
						{
							// physical_edges[abs(ee)]++;
						}
					}			
					break;
				}
				case 1:
				{
					auto vol1= abs(*vols.begin());
					/*if (domains[vol1] == conductor_id)
					{
						std::cout << "Possibile?" << std::endl;
						// intersurface[k]++;
						for (auto ee : _fte_list[k])
						{
							// physical_edges[abs(ee)]++;
						}
					}*/
					break;
				}
				case 0:
				{
					throw std::invalid_argument("Conductor boundary cannot be on mesh boundary!");
					break;
				}
			}	
		}
		
		std::vector<uint32_t>().swap(e_labels);
		_nte_list.resize(pts.size());
		// _etn_list.reserve(edges.size());
		// physical_nodes.resize(pts.size(),0);
		
		for (uint32_t k = 0; k < edges.size(); k++)
		{
			auto t = edges[k];
			
			uint32_t       p0(std::get<0>(t));
			uint32_t       p1(std::get<1>(t));
			
			sgnint32_t<int32_t> n1(p0,-1);
			sgnint32_t<int32_t> n2(p1, 1);
			
			std::vector<sgnint32_t<int32_t>> dummy(2);		
			_etn_list.push_back(dummy);
			_etn_list[k][0] = n1;
			_etn_list[k][1] = n2;
			
			sgnint32_t<int32_t> e1(k,-1);
			sgnint32_t<int32_t> e2(k, 1);
			
			_nte_list[p0].push_back(e1);
			_nte_list[p1].push_back(e2);
			
			// if (physical_edges[k])
			// {
				// for (auto nn : _etn_list[k])
				// {
					// physical_nodes[abs(nn)]++;
				// }
			// }
		}
		
		tc.toc();
		
		std::cout << "done - " << tc << " seconds" << std::endl;
		
		/************************ Read boundary surfaces ************************/
		linecount = 0;
		auto num_of_tets=lines;
		lines = strtot<uint32_t>(endptr, &endptr);
		
		tc.tic();
		while (linecount < lines)
		{
			/*if (ifs.fail())
			{
				std::cout << "Error while reading boundary surfaces" << std::endl;
				return false;
			}*/
			
			if ( (linecount%50000) == 0 )
			{
				std::cout << "Reading triangle: " << linecount;
				std::cout << "/" << lines << "\r";
				std::cout.flush();
			}
			
			auto t = parser::read_triangle_line<uint32_t>(endptr, &endptr);

			
			uint32_t       p0( std::get<1>(t) );
			uint32_t       p1( std::get<2>(t) );
			uint32_t       p2( std::get<3>(t) );
			uint32_t       bid( std::get<0>(t) );
			
			surface_type   tri( p0, p1, p2 );
			
			// if (!physical_surfaces[bid].size())
				// physical_surfaces[bid].resize(surfaces.size(),0);

			// auto itor = std::lower_bound(surfaces.begin(),surfaces.end(),tri);
			// physical_surfaces[bid][std::distance(surfaces.begin(),itor)]++;
			
			linecount++;
		}
		tc.toc();
		
		std::cout << "Reading triangle: " << linecount;
		std::cout << "/" << lines  << " - " << tc << " seconds"  << std::endl;

		tctot.toc();
		// std::cout << cyan << "Total time spent in reading mesh: ";
		// std::cout << tctot << " seconds" << nocolor << std::endl;
		
		// for (const auto& nn : _vtf_list)
			// assert(nn.size() == 4);
		
		// for (const auto& nn : _ftv_list)
			// assert(nn.size() == 2 || nn.size() == 1);
		
		// for (const auto& nn : _fte_list)
			// assert(nn.size() == 3);
		
		// for (const auto& nn : _etf_list)
		// {
			// assert(nn.size() >= 2);
		// }
		
		// for (const auto& nn : _nte_list)
			// assert(nn.size() >= 3);
		// for (const auto& nn : _etn_list)
			// assert(nn.size() == 2);	
		
		return true;
	}

	const std::vector<sgnint32_t<int32_t>>& vtf(const int32_t& v_id) const
	{
		if ( _vtf_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return _vtf_list[uint32_t(abs(v_id))];
	}

	// bool is_conductor(const uint32_t& v)
	// {
		// return domains[v] == conductor_id? true : false;
	// }

	const std::vector<sgnint32_t<int32_t>>& fte(const int32_t& f_id) const
	{
		if ( _fte_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return _fte_list[uint32_t(abs(f_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& ftv(const int32_t& f_id) const
	{
		if ( _ftv_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return _ftv_list[uint32_t(abs(f_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& etf(const int32_t& e_id) const
	{
		if ( _etf_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return _etf_list[uint32_t(abs(e_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& etn(const int32_t& e_id) const
	{
		if ( _etn_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return _etn_list[uint32_t(abs(e_id))];
	}

	const std::vector<sgnint32_t<int32_t>>& nte(const int32_t& n_id) const
	{
		if ( _nte_list.size() == 0 )
			throw std::invalid_argument("you have to populate the list first");
		else
			return _nte_list[uint32_t(abs(n_id))];
	}	

	Eigen::Vector3d face_barycenter(const uint32_t& f)
	{
		std::vector<uint32_t> nodes;
		Eigen::Vector3d bc(0);
		
		for (const auto& signed_ee : _fte_list[f])
		{
			uint32_t ee = signed_ee.Val();
			
			for (const auto& signed_nn : _etn_list[ee])
			{
				uint32_t nn = signed_nn.Val();
				
				if (!std::binary_search(nodes.begin(),nodes.end(),nn))
				{
					nodes.push_back(nn);
					bc += pts[nn];
				}	
			}
			
		}
		
		return bc/3;
	}

	Eigen::Vector3d vol_barycenter(const uint32_t& v)
	{
		std::vector<uint32_t> nodes;
		Eigen::Vector3d bc(0);
		
		for (const auto& signed_ff : _vtf_list[v])
		{
			uint32_t ff = signed_ff.Val();
			
			if (ff<surfaces.size())
			{
				for (const auto& signed_ee : _fte_list[ff])
				{
					uint32_t ee = signed_ee.Val();
					
					for (const auto& signed_nn : _etn_list[ee])
					{
						uint32_t nn = signed_nn.Val();
						
						if (!std::binary_search(nodes.begin(),nodes.end(),nn))
						{
							nodes.push_back(nn);
							bc += pts[nn];
						}	
					}
					
				}
			}
		}
		
		return bc/4;
	}

	Eigen::Vector3d edge_barycenter(const uint32_t& e)
	{
		Eigen::Vector3d bc(0);
		
		for (const auto& signed_nn : _etn_list[e])
		{
			uint32_t nn = signed_nn.Val();
			bc += pts[nn];	
		}
		
		return bc/2;
	}

	std::string print_face(const uint32_t& label, const uint32_t& f, bool orient, double r, double g, double b)
	{
		std::ostringstream fr;
		std::set<uint32_t> nodes;
		std::vector<Eigen::Vector3d > n;
		
		for (auto ee : fte(f))
		{
			for (auto nn : etn(abs(ee)))
				nodes.insert(abs(nn));
			
			if (nodes.size() == 3)
				break;
		}
		
		fr << "103.000 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		for (auto nn : nodes)
			n.push_back(pts[nn]);
		
		if (!orient)
			std::swap(n[1],n[2]);
		
		fr << n[0][0] << " " << n[0][1] << " " << n[0][2] << " " ;
		fr << n[1][0] << " " << n[1][1] << " " << n[1][2] << " " ;
		fr << n[2][0] << " " << n[2][1] << " " << n[2][2] << " " ;
			
		fr << std::endl;
		
		return fr.str();
	}

	std::string print_dual_edge(const uint32_t& label, const uint32_t& v, const uint32_t& f, bool orient, double r, double g, double b)
	{
		std::ostringstream fr;
		
		auto fb =  face_barycenter(f);
		auto vb =  vol_barycenter(v);
		
		fr << "102.100 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		if (orient)
			fr << fb[0] << " " << fb[1] << " " << fb[2] << " " << vb[0] << " " << vb[1] << " " << vb[2] << std::endl;
		else
			fr << vb[0] << " " << vb[1] << " " << vb[2] << " " << fb[0] << " " << fb[1] << " " << fb[2] << std::endl;
		
		return fr.str();
	}

	std::vector<std::string> print_dual_face(const uint32_t& label, const uint32_t& e, bool orient, double r, double g, double b)
	{
		// std::ostringstream fr;
		std::vector<std::string> ret;
		auto eb =  edge_barycenter(e);
		double x,y,z;
		x=y=z=0;
		int32_t new_or = orient? 1 : -1;
		
		
		for (auto n : _etn_list[e])
		{
			auto segno = n.Sgn();
			x+= new_or*segno*pts[abs(n)][0];
			y+= new_or*segno*pts[abs(n)][1];
			z+= new_or*segno*pts[abs(n)][2];
		}

		Eigen::Vector3d v1 { x, y, z };
		
		for (auto f : _etf_list[e])
		{
			auto fb = face_barycenter(abs(f));
			
			for (auto v : _ftv_list[abs(f)])
			{
				orient = new_or*f.Sgn()>0;
				auto vb =  vol_barycenter(abs(v));

				Eigen::Vector3d v2 = eb-vb;
				Eigen::Vector3d v3 = fb-vb;
				
				if (v1.dot(v2.cross(v3))>0)
					orient= true;
				else
					orient = false;
					
				
				ret.push_back(print_face(label,orient,vb[0],vb[1],vb[2],eb[0],eb[1],eb[2],fb[0],fb[1],fb[2],r,g,b));
			}
		}
		
		return ret;
	}

	std::string print_face(const uint32_t& f, const int32_t& orient)
	{
		std::ostringstream fr;
		std::set<uint32_t> nodes;
		
		for (auto ee : fte(f))
		{
			for (auto nn : etn(abs(ee)))
				nodes.insert(abs(nn));
			
			if (nodes.size() == 3)
				break;
		}
		
		for (auto nn : nodes)
			fr << nn << " " ;
			
		fr << orient << std::endl;
		
		return fr.str();
	}

	std::string print_edge(const uint32_t& label, const uint32_t& e, bool orient, double r, double g, double b)
	{
		std::ostringstream fr;
		
		const auto& nn = _etn_list[e];
		uint32_t nn_b = abs(*nn.begin());
		uint32_t nn_e = abs(*std::prev(nn.end()));
		Eigen::Vector3d n1 =  pts[nn_b];
		Eigen::Vector3d n2 =  pts[nn_e];
		
		fr << "102.100 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		if (orient)
			fr << n1[0] << " " << n1[1] << " " << n1[2] << " " << n2[0] << " " << n2[1] << " " << n2[2] << std::endl;
		else
			fr << n2[0] << " " << n2[1] << " " << n2[2] << " " << n1[0] << " " << n1[1] << " " << n1[2] << std::endl;
		
		return fr.str();
	}

	std::string print_edge(const uint32_t label, double x1, double y1, double z1, double x2, double y2, double z2)
	{
		std::ostringstream fr;
		
		fr << "102.000 " << label << " " << 0 << " " << 0 << " " << 0 << " 0.0 ";
		fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << std::endl;
		
		return fr.str();
	}

	std::string print_face(const uint32_t label, bool orient, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double r, double g, double b)
	{
		std::ostringstream fr;
		
		fr << "103.000 " << label << " " << r << " " << g << " " << b << " 0.0 ";
		if (orient)
			fr << x1 << " " << y1 << " " << z1 << " " << x2 << " " << y2 << " " << z2 << " " << x3 << " " << y3 << " " << z3 << std::endl;
		else
			fr << x1 << " " << y1 << " " << z1 << " " << x3 << " " << y3 << " " << z3 << " " << x2 << " " << y2 << " " << z2 << std::endl;
		
		return fr.str();
	}

	size_t volumes_size()
	{ 
		return volumes.size(); 
	}
	size_t surfaces_size()
	{ 
		return surfaces.size(); 
	}
	size_t edges_size()
	{ 
		return edges.size(); 
	}
	size_t nodes_size()
	{ 
		return pts.size(); 
	}
	
	private:
	std::vector<uint32_t> 						domains;
	std::vector<volume_type> 					volumes;
	std::vector<surface_type> 					surfaces;
	std::vector<edge_type> 						edges;
	std::vector<Eigen::Vector3d>	 			pts;
	/* node -> cluster of edge IDs around it */
	std::vector<cluster_list>    				_nte_list;
	std::vector<cluster_list>    				_etn_list;	
	/* edge -> cluster of face IDs around it */
	std::vector<cluster_list>        			_etf_list;
	std::vector<cluster_list>        			_fte_list;	
	/* triangle -> cluster of volume IDs around it */
	std::vector<cluster_list>  					_ftv_list;
	std::vector<cluster_list> 					_vtf_list;
};
