#include "FindInteractions.h"
//#include "../format.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>


using std::cout;
using std::endl;
using std::cerr;


void Binary::FindInteractions::readInputFile(const std::filesystem::path& path)
{
	std::ifstream inpFile(path.string());
	std::string line;
	if(inpFile.fail())
	{
		if(!path.parent_path().empty())
		{
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), path.relative_path().string()) << endl;
			cout << "Please check the file '"<< path.filename().string() <<"' in the followng directory\n";
			cout << "'" << path.relative_path().string() << "'" << endl;
		}
		else
		{
			std::filesystem::path p = std::filesystem::current_path();
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), p.relative_path().string()) << endl;
			cout << "Please check the file '"<< path.filename().string() <<"' in the followng directory\n";
			cout << "'" << p.relative_path().string() << "'" << endl;
		}
		exit(EXIT_FAILURE);
	}
	
	
	std::getline(inpFile, line);
	std::istringstream iss(line);
	std::string chain;
	while(iss >> chain)
	{
		binaryInterface += chain;
		m_vChains_binary.emplace_back(chain);
	}
	
}
