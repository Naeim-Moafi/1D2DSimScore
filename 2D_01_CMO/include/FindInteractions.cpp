#include "FindInteractions.h"
//#include "../format.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <chrono>

//using namespace std;
using namespace std::chrono;

using std::cout;
using std::endl;
using std::cerr;

void initMatrix(Matrix& mat, size_t length)
{
	mat.clear();
	for(size_t i { 0 }; i < length; ++i)
	{
		std::vector<char> v_tmp;
		for(size_t j { 0 }; j < length; ++j)
		{
			v_tmp.push_back('0');
		}
		mat.emplace_back(v_tmp);
	}
}

void CMO::FindInteractions::readInputFile(const std::filesystem::path& path)
{
	std::ifstream inpFile(path.string());
	std::string line;
	std::vector<std::string> v_lines;
	
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
	while(!inpFile.fail())
	{
		std::getline(inpFile, line);
		if(line.size() == 0)
		{ // break for emty line
			break;
		}
		v_lines.emplace_back(line);
	}
	
	initMatrix(m_Matrix,v_lines.size());
	//for(const auto& l : v_lines)
	//{
		//cout << l.size() << endl;
	//}
	//getchar();
	for(size_t i {0}; i < v_lines.size(); ++i)
	{
		for(size_t j {0}; j < v_lines.size(); ++j)
		{
			m_Matrix[i][j] = v_lines[i][j];
		}
	}
}


//#define __MAIN__
#ifdef __MAIN__
using namespace std;
//using namespace SS;
using namespace filesystem;
int main()
{
	CMO::FindInteractions fi;
	fi.readInputFile("3wbm.map");
}

#endif //__MAIN__
