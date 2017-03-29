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
 
 //file Mesh.hpp
 #include "Utilities.hpp"
 
class Mesh
{
	public:
	Mesh() 
	: type("none"), mesher("none")
	{
		loaded = false;
		xstep = ystep = zstep = 0;
	}

	void SetParam(uint32_t input_line, std::string param, std::string value)
	{
		if (param == "file")
			file = value;
		else if (param == "name")
			name = value;
		else if (param == "type")
		{
			if (std::find(meshtypes.begin(),meshtypes.end(),value) == meshtypes.end())
				MyThrow(input_line,mesh_unknown_type);
			type = value;
		}
		else if (param == "mesher")
		{
			if (std::find(meshers.begin(),meshers.end(),value) == meshers.end())
				MyThrow(input_line,mesh_unknown_mesher);
			mesher = value;
		}
		else if (param == "xstep")
			xstep = std::stod(value);
		else if (param == "ystep")
			ystep = std::stod(value);
		else if (param == "zstep")
			zstep = std::stod(value);
		else
			MyThrow(input_line,mesh_unknown_parameter);
	}
	
	const std::string& Name() { return name; }
	const std::string& FileName() { return file; }
	const std::string& GetMeshType() { return type; }
	const std::string& GetMesher()   { return mesher; }
	const double& GetLx() { return xstep; }
	const double& GetLy() { return xstep; }
	const double& GetLz() { return xstep; }
	bool IsLoaded() { return loaded; }
	void Switch() { loaded = !loaded; }
	
	private:
	std::string file, name, type;
	std::string mesher;
	bool loaded;
	double xstep,ystep,zstep; //used only when mesh type is cartesian
};
