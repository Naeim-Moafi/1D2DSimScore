#include <algorithm>
#include <exception>
#include "FindInteractions.h"
//#include "../format.h"
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::cerr;


SS_all::FindInteractions::FindInteractions()
	: m_numberOfCanBasePairs{0}, m_numberOfNonCanBasePairs{0}, m_numberOfAllBasePairs{0} 
	{
		m_isWobble_canonical = false;
	}

SS_all::FindInteractions::FindInteractions(bool isWobble_canonical)
	: FindInteractions()
{
	m_isWobble_canonical = isWobble_canonical;
}

bool SS_all::FindInteractions::get_isWobble_canonical() const
{
	return m_isWobble_canonical;
}

void SS_all::FindInteractions::set_isWobble_canonical(bool isWobble_canonical)
{
	SS::FindInteractions::set_isWobble_canonical(isWobble_canonical);
}

std::string SS_all::FindInteractions::get_extension() const
{
	return m_extension;
}

void SS_all::FindInteractions::set_extension(const std::string& extension)
{
	m_extension = extension;
}


void SS_all::FindInteractions::init_sequence(const std::filesystem::path& path)
{
	if(std::filesystem::is_regular_file(path))
	{
		//cout << "read from file directly\n";
		// it means you can find everything in a single files
		std::ifstream inpFile(path.string());
		std::string line;
		// check the file
		if(inpFile.fail())
		{
			if(!path.parent_path().empty())
			{
				//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), path.relative_path().string()) << endl;
				cout << "please check the file '" << path.filename().string() << "' in the following directory\n";
				cout << "'" << path.relative_path().string() << "'" << endl;
			}
			else
			{
				std::filesystem::path p = std::filesystem::current_path();
				//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), p.relative_path().string()) << endl;
				cout << "please check the file '" << path.filename().string() << "' in the following directory\n";
				cout << "'" << p.relative_path().string() << "'" << endl;
			}
			exit(EXIT_FAILURE);
		}
		
		while(std::getline(inpFile, line))
		{
			if(line == ">seq")
			{
				std::getline(inpFile, line);
				m_sequence = line;
				SS::FindInteractions::m_sequence = m_sequence;
				m_seqFlag = true;
				break;
			}
		}
		inpFile.close();
	}
	else if(std::filesystem::is_directory(path))
	{
		//cout << "finding sequence file in folder and read it\n";
		const std::filesystem::path target_dir { path };
		// iterate all over the target directory to find the sequence file
		for(const auto& dir_entry : std::filesystem::directory_iterator { target_dir })
		{
			const std::filesystem::path tmp_path { dir_entry };
			if(tmp_path.extension() == ".seq")
			{
				SS::FindInteractions::init_sequence(tmp_path);
				std::ostringstream ossSeq;
				m_sequence = SS::FindInteractions::m_sequence;
				m_seqFlag = true;
			}
		}
	}
	else
	{
		cout << "12DDSimScore does not support this type of the input\n";
		cout << "use a regular file or a folder as input" << endl;
		exit(EXIT_FAILURE);
	}
	
	//std::cout << m_sequence << std::endl;
	//getchar();
}

void SS_all::FindInteractions::init_structures_single(const std::filesystem::path& path)
{
	std::ifstream inpFile(path.string());
	std::string line;
	SSMap ssm;
	SSMatrix ssMat;
	init_sequence(path);
	// check if file is opened successfully
	if(inpFile.fail())
	{
		if(!path.parent_path().empty())
		{
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), path.relative_path().string()) << endl;
			cout << "please check the file '" << path.filename().string() << "' in the following directory\n";
			cout << "'" << path.relative_path().string() << "'" << endl;
		}
		else
		{
			std::filesystem::path p = std::filesystem::current_path();
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), p.relative_path().string()) << endl;
			cout << "please check the file '" << path.filename().string() << "' in the following directory\n";
			cout << "'" << p.relative_path().string() << "'" << endl;
		}
		exit(EXIT_FAILURE);
	}
	
	while(std::getline(inpFile, line))
	{
		if(line.size() == 0)
		{
			break;
		}
		cout << line << endl;
		if(line[0] == '>' && line != ">seq")
		{
			m_vStructures_name.push_back(line.substr(1));
			auto oldPos = inpFile.tellg();
			std::vector<std::string> v_lines;
			while(std::getline(inpFile, line))
			{
				if(line.size() == 0)
				{
					break;
				}
				if(line[0] == '>')
				{
					inpFile.seekg(oldPos);
					break;
				}
				//cout << line << endl;
				v_lines.push_back(line);
			}
			m_vSSMaps.emplace_back(SS::FindInteractions::extractInteractions_vector(v_lines));
			m_vMatrix_all.emplace_back(SS::FindInteractions::extractInteractions_matrix(v_lines));
		}
	}
}

void SS_all::FindInteractions::init_structures_multiple(const std::filesystem::path& path)
{
	const std::filesystem::path target_dir { path };
	init_sequence(path);
	for(const auto& dir_entry : std::filesystem::directory_iterator { target_dir })
	{
		const std::filesystem::path tmp_path { dir_entry };
		if(tmp_path.extension() == m_extension)
		{
			m_extFlag = true;
			std::ifstream inpFile(tmp_path.string());
			m_vStructures_name.emplace_back(tmp_path.stem().string());
			std::string line;
			std::vector<std::string> v_lines;
			while(std::getline(inpFile, line))
			{
				if(line.size() == 0)
				{
					break;
				}
				v_lines.push_back(line);
			}
			m_vSSMaps.emplace_back(SS::FindInteractions::extractInteractions_vector(v_lines));
			m_vMatrix_all.emplace_back(SS::FindInteractions::extractInteractions_matrix(v_lines));
		}
	}
	
	if(!m_extFlag)
	{
		cout << "There is no file with " << m_extension << " in your dataset\n";
		exit(EXIT_FAILURE);
	}
}

void SS_all::FindInteractions::init_structures(const std::filesystem::path& path)
{
	
	std::vector<std::vector<int>> matrix_map;
	SSContacts ssc;
	SSMap ssm;
	SSMatrix ssMat;
	// check if path is a regular file
	if(std::filesystem::is_regular_file(path))
	{
		init_structures_single(path);
	}
	else if(std::filesystem::is_directory(path))
	{
		init_structures_multiple(path);
	}
	else
	{
		cout << "1DSimScore does not support this type of the input\n";
		cout << "use a regular file or a folder as input" << endl;
		exit(EXIT_FAILURE);
	}
	
	
}




//#define __MAIN__
#ifdef __MAIN__
using namespace std;
//using namespace SS;
using namespace filesystem;
int main()
{	
	SS_all::FindInteractions findSSAll;
	//findSSAll.init_sequence("freeSL2");
	//findSSAll.init_sequence("samples/freeSL2/AllInOne.SS_all");
	//findSSAll.init_structures("freeSL2/AllInOne.SS_all");
	findSSAll.init_structures("samples/freeSL2/AllInOne.SS_all");
	return 0;
}

#endif //__MAIN__
