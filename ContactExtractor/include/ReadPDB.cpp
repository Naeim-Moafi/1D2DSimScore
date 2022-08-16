#include <algorithm>
#include "ReadPDB.h"
#include <iomanip>
#include <sstream>

const int MAX_N_ITEMS = 4;
const int MAX_STR_LENGTH = 32;

bool are_there_spaces(const std::string& inp_str)
{
	auto itSpace =  std::find(inp_str.begin(), inp_str.end(), ' ');
	if(itSpace != inp_str.end())
	{
		return true;
	}
	
	return false;
}

std::string del_spaces(const std::string& inp_str)
{
	std::string out_str;
	for(size_t i = 0; i < inp_str.size(); ++i)
	{
		if(inp_str[i] == ' ')
		{
			continue;
		}
		out_str += inp_str[i];
	}
	
	return out_str;
}

int split_string(std::vector<std::string>& output_vector, const std::string& string2parse, char delim)
{
	output_vector.clear();
	std::istringstream iss(string2parse);
	std::string curr_token;
	int n_items = 0;
	while(std::getline(iss, curr_token, delim))
	{ 
		++n_items;
		output_vector.emplace_back(curr_token);
	}
	
	//if(n_items >= MAX_N_ITEMS)
	//{
		//cout << "too many items\n";
		//exit(EXIT_FAILURE);
	//}
	
	return n_items;
}


void ReadPDB::read_file(const std::filesystem::path& pdbPath)
{
	std::ifstream inpFile(pdbPath.string());
	std::string line;
	n_atoms = 0;
	n_residues = 0;
	n_chains = 0;
	Atom curr_atom;
	std::vector<Atom> curr_chain;
	if(inpFile.fail())
	{
		if(!pdbPath.parent_path().empty())
		{
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), path.relative_path().string()) << endl;
			cout << "Please check the file '" << pdbPath.filename().string() << "' in the followng directory\n'" << pdbPath.relative_path().string() << "'" << endl;
		}
		else
		{
			std::filesystem::path p = std::filesystem::current_path();
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), p.relative_path().string()) << endl;
			cout << "Please check the file '" << p.filename().string() << "' in the followng directory\n'" << p.relative_path().string() << "'" << endl;
		}
		exit(EXIT_FAILURE);
	}
	while(std::getline(inpFile, line))
	{
		if(line.substr(0,3) == "TER")
		{
			++n_chains;
			m_vAtomsInChains.emplace_back(curr_chain);
			curr_chain.clear();
			continue;
		}
		if(line == "ENDMDL" || del_spaces(line) == "END" || line.size() == 0)
		{
			if(curr_chain.size() > 0)
			{
				++n_chains;
				m_vAtomsInChains.emplace_back(curr_chain);
				curr_chain.clear();
			}
			break;
		}
		
		if(fill_atom_strc(curr_atom, line))
		{
			++n_atoms;
			m_vAtoms.emplace_back(curr_atom);
			curr_chain.emplace_back(curr_atom);
		}
	}
	//m_vResidues.reserve(n_atoms);
	//m_vChains.resize(n_chains);
	//m_vResidues[0].start_index = 0;
	
	//int n_atoms_in_res = 0;
	
	//for(int i{0}; i < n_atoms-1; ++i)
	//{
		//++n_atoms_in_res;
		//if(m_vAtoms[i].nuclNumber != m_vAtoms[i+1].nuclNumber)
		//{
			//m_vResidues[n_residues].n_atoms = n_atoms_in_res;
			//m_vResidues[n_residues+1].start_index = i + 1;
			//++n_residues;
			//n_atoms_in_res = 0;
		//}
	//}
	//m_vResidues[n_residues].n_atoms = n_atoms_in_res+1;
	//++n_residues;
	//int n_residues_in_chain = 0;
	//m_vChains[0].start_index = 0;
	//m_vChains[0].n_residues_sum = 0;
	//int residues_counter = 0;
	//int i_chain = 0;
	//for(int i{0}; i < n_residues-1; ++i)
	//{
		//++n_residues_in_chain;
		//++residues_counter;
		//if(m_vAtoms[m_vResidues[i].start_index].chainId != m_vAtoms[m_vResidues[i+1].start_index].chainId)
		//{
			//m_vChains[i_chain].n_residues = n_residues_in_chain;
			//m_vChains[i_chain+1].n_residues_sum = residues_counter;
			//m_vChains[i_chain+1].start_index = m_vResidues[i+1].start_index;
			//++i_chain;
			//n_residues_in_chain = 0;
		//}
	//}
	//m_vChains[i_chain].n_residues = n_residues_in_chain + 1;
	//selected_chain=0;
}

bool ReadPDB::fill_atom_strc(Atom& atom, const std::string& line)
{
	std::string temp_buffer;
	[[maybe_unused]]int inputLength = line.size();
	
	temp_buffer = line.substr(0,6);
	if(temp_buffer != "ATOM  " && temp_buffer != "HETATM") 
	{
		return false;
	}
	
	try
	{
		temp_buffer = line.substr(6,5);
	}
	catch(...)
	{
		return false;
	}
	
	atom.atomNumber = std::stoi(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(12, 4);
	}
	catch(...)
	{
		return false;
	}
	
	atom.atomName = temp_buffer;
	
	try
	{
		temp_buffer = line.substr(17,3);
	}
	catch(...)
	{
		return false;
	}
	
	atom.nuclName = temp_buffer;
	
	try
	{
		temp_buffer = line.substr(20,2);
	}
	catch(...)
	{
		return false;
	}
	
	atom.chainId = temp_buffer;
	
	try
	{
		temp_buffer = line.substr(22,4);
	}
	catch(...)
	{
		return false;
	}
	
	atom.nuclNumber = std::stoi(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(30,8);
	}
	catch(...)
	{
		return false;
	}
	
	atom.coord.x = std::stod(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(38,8);
	}
	catch(...)
	{
		return false;
	}
	
	atom.coord.y = std::stod(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(46,8);
	}
	catch(...)
	{
		return false;
	}
	
	atom.coord.z = std::stod(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(54,6);
	}
	catch(...)
	{
		temp_buffer = "1.0";
	}
	
	atom.occupancy = std::stod(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(54,6);
	}
	catch(...)
	{
		temp_buffer = "1.0";
	}
	
	atom.occupancy = std::stod(temp_buffer);
	
	try
	{
		temp_buffer = line.substr(60,6);
	}
	catch(...)
	{
		temp_buffer = "0.0";
	}
	
	atom.bfactor = std::stod(temp_buffer);
	
	return true;
}

int ReadPDB::get_n_atoms() const
{
	return n_atoms;
}
