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
		if(line == "ENDMDL" || del_spaces(line) == "END")
		{
			break;
		}
		++n_atoms;
		if(fill_atom_strc(curr_atom, line))
		{
			m_vAtoms.emplace_back(curr_atom);
			curr_chain.emplace_back(curr_atom);
		}
	}
	
	m_vResidues.reserve(n_atoms);
	m_vChains.resize(n_chains);
	m_vResidues[0].start_index = 0;
	
	
	int n_atoms_in_res = 0;
	
	for(int i{0}; i < n_atoms-1; ++i)
	{
		++n_atoms_in_res;
		if(m_vAtoms[i].nuclNumber != m_vAtoms[i+1].nuclNumber)
		{
			m_vResidues[n_residues].n_atoms = n_atoms_in_res;
			m_vResidues[n_residues+1].start_index = i + 1;
			++n_residues;
			n_atoms_in_res = 0;
		}
	}
	
	m_vResidues[n_residues].n_atoms = n_atoms_in_res+1;
	++n_residues;
	
	int n_residues_in_chain = 0;
	m_vChains[0].start_index = 0;
	m_vChains[0].n_residues_sum = 0;
	int residues_counter = 0;
	int i_chain = 0;
	
	for(int i{0}; i < n_residues-1; ++i)
	{
		++n_residues_in_chain;
		++residues_counter;
		if(m_vAtoms[m_vResidues[i].start_index].chainId != m_vAtoms[m_vResidues[i+1].start_index].chainId)
		{
			m_vChains[i_chain].n_residues = n_residues_in_chain;
			m_vChains[i_chain+1].n_residues_sum = residues_counter;
			m_vChains[i_chain+1].start_index = m_vResidues[i+1].start_index;
			++i_chain;
			n_residues_in_chain = 0;
		}
	}
	
	m_vChains[i_chain].n_residues = n_residues_in_chain + 1;
	

	//for(const auto& curr_chain_atoms : m_vAtomsInChains)
	//{
		//int n_residues_in_chain = 0;
		//for(size_t i = 0; i < curr_chain_atoms.size(); ++i)
		//{
			//++n_atoms_in_res;
			//if(curr_chain_atoms[i].nuclNumber != curr_chain_atoms[i+1].nuclNumber)
			//{
				//m_vResidues[n_residues].n_atoms = n_atoms_in_res;
				//m_vResidues[n_residues+1].start_index = i + 1;
				//cout << n_residues" " << m_vResidues[n_residues].start_index << " " << m_vResidues[n_residues+1].start_index << endl;
				//++n_residues;
				//++n_residues_in_chain;	
				//n_atoms_in_res = 0;
			//}
		//}
		//if(i_chain == 0)
		//{
			//m_vChains[i_chain].n_residues = n_residues_in_chain;
			//m_vChains[i_chain].n_residues_sum = n_residues;
			//m_vChains[i_chain].start_index = n_residues - n_residues_in_chain;
			//++i_chain;
			//continue;
		//}
		//m_vChains[i_chain].n_residues = n_residues_in_chain;
		//m_vChains[i_chain].n_residues_sum = n_residues;
		//m_vChains[i_chain].start_index = n_residues - n_residues_in_chain;
		//++i_chain;
	//}
	selected_chain=0;
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

int ReadPDB::get_n_chains()const
{
	return n_chains;
}

int ReadPDB::get_n_residues()const
{
	return n_residues;
}

int ReadPDB::select_chain_for_reading(int i_chain)
{
	if(i_chain>=n_chains || i_chain<0)
	{
		cout << "ReadPDB::" << __func__ << ": reference to the chain " <<  i_chain  << " which is outside the table size 0-" << n_chains-1 << endl;
		exit(EXIT_FAILURE);
	}
	
	
	selected_chain = i_chain;
	return selected_chain;
}

int ReadPDB::get_n_residues_in_chain(int i_chain) const
{
	if(i_chain>=n_chains || i_chain<0)
	{
		cout << "ReadPDB::" << __func__ << ": reference to the chain " <<  i_chain  << " which is outside the table size 0-" << n_chains-1 << endl;
		exit(EXIT_FAILURE);
	}
	
	return m_vChains[i_chain].n_residues;
}

int ReadPDB::get_atom(Atom& i_atom, int i_res, const std::string& i_name)const
{
	int n_atoms_found, i, j, start_index, n_atoms_in_res, n_splitted_elems;
    std::string tmp_str_name, tmp_str_to_split;
    std::vector<std::string> str_splitted;
    str_splitted.reserve(MAX_N_ITEMS);
    
    if(i_res>=n_residues || i_res<0)
    {
		cout << "ReadPDB::" << __func__ << ": request getting residue " << i_res << " which is out of the residue table size 0-" << n_residues-1 << endl;
		exit(EXIT_FAILURE);
    }
	start_index = m_vResidues[i_res].start_index;
	n_atoms_in_res = m_vResidues[i_res].n_atoms;
	n_atoms_found = 0;
	for( i = start_index; i < start_index + n_atoms_in_res; ++i)
	{
		if(are_there_spaces(i_name))
		{
			tmp_str_name = m_vAtoms[i].atomName;
			tmp_str_to_split = i_name.substr(0,MAX_STR_LENGTH);
			str_splitted[0] = tmp_str_to_split;
			n_splitted_elems = 1;
		}
		else
		{
			tmp_str_name = del_spaces(m_vAtoms[i].atomName);
			tmp_str_to_split = i_name.substr(0,MAX_STR_LENGTH);
			n_splitted_elems = split_string(str_splitted, tmp_str_to_split,'|');
		}
		for(j = 0; j  < n_splitted_elems; ++j)
		{
			if(str_splitted[j] == tmp_str_name)
			{
				if(n_atoms_found == 0)
				{
					i_atom.atomNumber=m_vAtoms[i].atomNumber;
					i_atom.nuclNumber=m_vAtoms[i].nuclNumber;
					i_atom.atomName = del_spaces(m_vAtoms[i].atomName);
					i_atom.nuclName = del_spaces(m_vAtoms[i].nuclName);
					i_atom.chainId= del_spaces(m_vAtoms[i].chainId);
					i_atom.coord.x=m_vAtoms[i].coord.x;
					i_atom.coord.y=m_vAtoms[i].coord.y;
					i_atom.coord.z=m_vAtoms[i].coord.z;
					
					i_atom.coordBackup.x=m_vAtoms[i].coord.x;
					i_atom.coordBackup.y=m_vAtoms[i].coord.y;
					i_atom.coordBackup.z=m_vAtoms[i].coord.z;
					
					i_atom.occupancy = m_vAtoms[i].occupancy;
					i_atom.bfactor = m_vAtoms[i].bfactor;

					n_atoms_found++;
				}
				else
				{
					++n_atoms_found;
				}
			}
		}
	}
	
	if(n_atoms_found != 1)
	{

		if(n_atoms_found == 0)
		{
			if(n_atoms_in_res > 0)
			{
				i = start_index;
			}
			else
			{
				i_atom.atomName.clear();
			}
		}
		else
		{
			i = start_index;
			cout << "more than 1 atom (" << n_atoms_found << " atoms found) _" << i_name << "_ in residue " << m_vAtoms[i].nuclNumber << " chain " << m_vAtoms[i].chainId << "\n";
		}
	}
	return n_atoms_found;
}

int ReadPDB::get_atom(Atom& i_atom, int i_res, const std::string& i_name, int i_chain)const
{
	if(i_chain>=n_chains || i_chain<0)
	{
		cout << "ReadPDB::" << __func__ << ": reference to the chain " <<  i_chain  << " which is outside the table size 0-" << n_chains-1 << endl;
		exit(EXIT_FAILURE);
	}
	
    if(i_res>=n_residues || i_res<0)
    {
		cout << "ReadPDB::" << __func__ << ": request getting residue " << i_res << " which is out of the residue table size 0-" 
		<< m_vChains[i_chain].n_residues-1 << " for chain: " << i_chain << endl;
		exit(EXIT_FAILURE);
    }
	
	
	return get_atom(i_atom, m_vChains[i_chain].n_residues_sum + i_res, i_name);
}

int ReadPDB::get_atom_from_selected_chain(Atom& i_atom, int i_res, const std::string& i_name) const
{
	if(i_res>=n_residues || i_res<0)
    {
		cout << "ReadPDB::" << __func__ << ": request getting residue " << i_res << " which is out of the residue table size 0-" 
		<< m_vChains[selected_chain].n_residues-1 << " for chain: " << selected_chain << endl;
		exit(EXIT_FAILURE);
    }
    
	return get_atom(i_atom, m_vChains[selected_chain].n_residues_sum + i_res, i_name);
}

//#define __MAIN__
#ifdef __MAIN__
using namespace std;
int main()
{
	ReadPDB inputPDB;
	inputPDB.read_file("sample.pdb");
	//cout << inputPDB.m_vChains.size() << endl;
	
	cout << inputPDB.get_n_atoms() << " " << inputPDB.get_n_residues() << " " << inputPDB.get_n_chains() << endl;
	//for(int i = 0; i < inputPDB.get_n_chains(); ++i)
	//{
		//cout << inputPDB.get_n_residues_in_chain(i+1) << endl;
	//}
	
	//string s("na|eim");
	int numberOfChains = inputPDB.get_n_chains();
	int numberOfNucleotides = inputPDB.get_n_residues();
	std::vector<Range> chains(numberOfChains);
	Atom nucleotides[numberOfNucleotides];
	std::vector<Atom> m_vNucleotides;
	m_vNucleotides.reserve(numberOfNucleotides);
	for(int i{0}; i < numberOfChains; ++i)
	{
		auto numberOfNucleotidesInCurrChain = inputPDB.get_n_residues_in_chain(i);
		inputPDB.select_chain_for_reading(i);
		for(int j{0}; j < numberOfNucleotidesInCurrChain; ++j)
		{
			// using atom C4' for extracting informations
			int isReadC4_OK = inputPDB.get_atom_from_selected_chain(nucleotides[chains[i].start + j], j, "C4'|C4*");
			// check if the pdb file has gaps or not
			m_vNucleotides.push_back(nucleotides[chains[i].start + j]);			
		}
	}
	int i = 0;
	for(auto const& atom : m_vNucleotides)
	{
		cout << std::setw(2) << ++i << " " << atom.nuclName  << " " << atom.nuclNumber << " " << atom.atomName  << " " << atom.atomNumber << atom.chainId << endl;
	}
	//for(auto const& atom : inputPDB.m_vAtoms)
	//{
		//cout << std::setw(2) << ++i << " " << atom.nuclName  << " " << atom.nuclNumber << " " << atom.atomName << atom.chainId << endl;
	//}
		
	//vector<string> v_ouput(MAX_N_ITEMS);
	//split_string(v_ouput, s, '|');
	//for(const auto& entry : v_ouput)
	//{
		//cout << entry << endl;
	//}
	
	
	
	//cout << s << endl;
	
	
	//getchar();();
	//for(const auto& atoms : inputPDB.m_vAtomsInChains[3])
	//{
		//cout << atoms.nuclNumber << endl;
	//}
	
	return 0;
}

#endif
