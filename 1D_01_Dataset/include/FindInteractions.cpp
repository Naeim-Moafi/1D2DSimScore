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


void::Binary_all::FindInteractions::set_extension(const std::string& extension)
{
	m_extension = extension;
}

std::string Binary_all::FindInteractions::get_extension() const
{
	return m_extension;
}

void Binary_all::FindInteractions::init_binary_single(const std::filesystem::path& path)
{
	std::ifstream inpFile(path.string());
	std::string line;
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
		if(line[0] == '>')
		{
			m_vNames.push_back(line.substr(1));
		}
		else
		{
			std::istringstream iss(line);
			std::string chain;
			std::string binaryInterface;
			while(iss >> chain)
			{
				binaryInterface += chain;
			}
			m_vBinaryInterface.emplace_back(binaryInterface);
		}
	}
}

void Binary_all::FindInteractions::init_binary_multiple(const std::filesystem::path& path)
{
	const std::filesystem::path target_dir { path };
	for(const auto& dir_entry : std::filesystem::directory_iterator { target_dir })
	{
		const std::filesystem::path tmp_path { dir_entry };
		if(tmp_path.extension() == m_extension)
		{
			m_extFlag = true;
			std::ifstream inpFile(tmp_path.string());
			m_vNames.emplace_back(tmp_path.stem().string());
			std::string line;
			std::getline(inpFile, line);
			std::istringstream iss(line);
			std::string chain;
			std::string binaryInterface;
			while(iss >> chain)
			{
				binaryInterface += chain;
			}
			m_vBinaryInterface.emplace_back(binaryInterface);
		}
	}
	
	if(!m_extFlag)
	{
		cout << "There is no file with " << m_extension << " in your dataset\n";
		exit(EXIT_FAILURE);
	}
}

void Binary_all::FindInteractions::init_binary(const std::filesystem::path& path)
{
	// check if path is a regular file
	if(std::filesystem::is_regular_file(path))
	{
		init_binary_single(path);
	}
	else if(std::filesystem::is_directory(path))
	{
		init_binary_multiple(path);
	}
	else
	{
		cout << "1DSimScore does not support this type of the input\n";
		cout << "use a regular file or a folder as input" << endl;
		exit(EXIT_FAILURE);
	}
}
