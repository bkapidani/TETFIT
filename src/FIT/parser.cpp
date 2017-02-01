//file parser.cpp

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <string>

int main(int argc, char **argv)
{
    assert(argc==2);
	
	std::ifstream preamble,postamble, geom;
	std::ofstream dest("mesher.hpp",std::ios::binary);
	
	preamble.open("preamble.cpp",std::ios::binary);
	dest << preamble.rdbuf();
	preamble.close();
	dest.close();
	
	dest.open("mesher.hpp",std::ofstream::out | std::ofstream::app);
	geom.open(argv[1]);
	
	std::string new_instr;
	std::string semicolon(";");
	std::vector<std::string> diseq_vec;
	for (size_t k=0; k<7; k++) //min and max values of x,y,z and grid resolution h
	{
		std::getline(geom, new_instr,'\n');
		auto instr = new_instr + semicolon;
		dest << "      " << instr << std::endl;
	}
	
	while (std::getline(geom,new_instr,'\n'))
		diseq_vec.push_back(new_instr);
	
	postamble.open("postamble.cpp",std::ios::binary);
	dest << postamble.rdbuf();
	postamble.close();
	
	for (auto d : diseq_vec)
	{
		dest << "      if (!( " << d << " ))" << std::endl;
		dest << "         return false;" << std::endl;
	}
	
	dest << "      return true;" << std::endl;
    dest << "   }" << std::endl << "};" << std::endl << "#endif" << std::endl;
	
	dest.close();
	geom.close();
	
	std::cout << "-- Parsing geometry from file done" << std::endl;
    return 0;
}
