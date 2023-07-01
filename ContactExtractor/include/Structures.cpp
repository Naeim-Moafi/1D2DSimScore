#include <fstream>
#include "Structures.h"
#include <iterator>
#include <sstream>
#include <thread>
#include <cstdio>
#include <iomanip>
std::vector<std::string> RNA = {"A", "C", "G", "U"};
std::vector<std::string> DNA = {"DA", "DC", "DG", "DU", "DT"};
std::vector<std::string> PROTEIN = {"ALA", "ARG", "ASN", "ASP",
									"CYS", "GLN", "GLU", "GLY",
									"HIS", "ILE", "LEU", "LYS",
									"MET", "PHE", "PRO", "SER",
									"THR", "TRP", "TYR", "VAL"};
									
std::map<std::string, std::string> DNA_SEQ = {{"DA", "A"}, {"DC", "C"}, {"DG", "G"}, {"DU", "U"}, {"DT", "T"}};	
									
std::map<std::string, std::string> PROTEIN_SEQ = {{"ALA", "A"}, {"ARG", "R"}, {"ASN", "N"}, {"ASP", "D"},
												  {"CYS", "C"}, {"GLN", "Q"}, {"GLU", "E"}, {"GLY", "G"},
												  {"HIS", "H"}, {"ILE", "I"}, {"LEU", "L"}, {"LYS", "K"}, 
												  {"MET", "M"}, {"PHE", "F"}, {"PRO", "P"}, {"SER", "S"}, 
												  {"THR", "T"}, {"TRP", "W"}, {"TYR", "Y"}, {"VAL", "V"}};
												  


void initial_conatact_map(std::vector<std::vector<int>>& contact_map, size_t size)
{
	for(size_t i { 0 }; i < size; ++i)
	{
		std::vector<int> v_tmp;
		for(size_t j { 0 }; j < size; ++j)
		{
			v_tmp.push_back(0);
		}
		contact_map.emplace_back(v_tmp);
	}
}

void initial_conatact_map(std::vector<std::vector<int>>& contact_map, size_t size1, size_t size2)
{
	for(size_t i { 0 }; i < size1; ++i)
	{
		std::vector<int> v_tmp;
		for(size_t j { 0 }; j < size2; ++j)
		{
			v_tmp.push_back(0);
		}
		contact_map.emplace_back(v_tmp);
	}
}


std::vector<std::string> structure_string2vec(const std::string& inp_struct, const std::vector<std::string>& v_seq)
{
	int n = 0;
	std::vector<std::string> v_struct;
	for(auto const& chain : v_seq)
	{
		v_struct.emplace_back(inp_struct.substr(n,chain.size()));
		n += chain.size();
	}
	
	return v_struct;
}

void clear_v_string(std::vector<std::string>& v_str)
{
	for(auto& str : v_str) { str.clear();}
}

void Structures::readPDB_allStructures(const std::filesystem::path& inputPath)
{
	ReadPDB inputPDB;
	inputPDB.read_file(inputPath);
	m_vAtoms_allStructures = inputPDB.m_vAtoms;
	sepRNA();
	sepDNA();
	sepProtein();
	//sepLigand();
	
	m_max_n_threads = 8;
	m_rna_seq = extract_rna_seq();
	m_dna_seq = extract_dna_seq();
	m_protein_seq = extract_protein_seq();
	create_all_seq();
	
	
	//update
	if(update)
	{
		read_restraint_file(selection_path);
		sepRNA();
		sepDNA();
		sepProtein();
		clear_v_string(m_rna_seq);
		clear_v_string(m_dna_seq);
		clear_v_string(m_protein_seq);
		clear_v_string(m_sequence_sep_chains);
		
		rna_size = dna_size = protein_size = 0;
		
		m_rna_seq = extract_rna_seq();
		m_dna_seq = extract_dna_seq();
		m_protein_seq = extract_protein_seq();
		create_all_seq();
		// add structures in desired order
		m_vAtoms_allStructures.clear();
		m_vAtoms_allStructures.insert(m_vAtoms_allStructures.end(), m_vAtoms_rna.begin(), m_vAtoms_rna.end());
		m_vAtoms_allStructures.insert(m_vAtoms_allStructures.end(), m_vAtoms_dna.begin(), m_vAtoms_dna.end());
		m_vAtoms_allStructures.insert(m_vAtoms_allStructures.end(), m_vAtoms_protein.begin(), m_vAtoms_protein.end());
	}
	else
	{
		// add structures in desired order
		m_vAtoms_allStructures.clear();
		m_vAtoms_allStructures.insert(m_vAtoms_allStructures.end(), m_vAtoms_rna.begin(), m_vAtoms_rna.end());
		m_vAtoms_allStructures.insert(m_vAtoms_allStructures.end(), m_vAtoms_dna.begin(), m_vAtoms_dna.end());
		m_vAtoms_allStructures.insert(m_vAtoms_allStructures.end(), m_vAtoms_protein.begin(), m_vAtoms_protein.end());
		

		std::vector<Atom> v_atoms;
		std::string chain_seq;
		Residue residue;
		Chain chain;
		for(size_t i_atom  = 0; i_atom < m_vAtoms_allStructures.size()-1; ++i_atom)
		{
			v_atoms.emplace_back(m_vAtoms_allStructures[i_atom]);
			if(m_vAtoms_allStructures[i_atom].nuclNumber != m_vAtoms_allStructures[i_atom + 1].nuclNumber)
			{
				residue.v_atoms_in_residue.insert(residue.v_atoms_in_residue.end(), v_atoms.begin(), v_atoms.end());
				v_atoms.clear();
				residue.residue_name = del_spaces(m_vAtoms_allStructures[i_atom].nuclName);
				chain.v_residues_in_chain.emplace_back(residue);
				residue.v_atoms_in_residue.clear();
				if(m_vAtoms_allStructures[i_atom].chainId != m_vAtoms_allStructures[i_atom + 1].chainId)
				{
					m_vChains.push_back(chain);
					chain.v_residues_in_chain.clear();
				}
			}
		}
		
		residue.v_atoms_in_residue.insert(residue.v_atoms_in_residue.end(), v_atoms.begin(), v_atoms.end());
		chain.v_residues_in_chain.emplace_back(residue);
		m_vChains.push_back(chain);
	}
	
	// initializing ranges;
	m_rna_range.start = m_dna_range.start = m_protein_range.start = 0;
	m_rna_range.end = m_dna_range.end = m_protein_range.end = 0;
	
	if(m_vAtoms_rna.size() > 0)
	{
		m_rna_range.start = 0;
		m_rna_range.end = rna_size - 1;
	}
	
	if(m_vAtoms_dna.size() > 0)
	{
		m_dna_range.start = rna_size;
		m_dna_range.end = m_dna_range.start + dna_size - 1;
	}
	
	if(m_vAtoms_protein.size() > 0)
	{
		m_protein_range.start = rna_size + dna_size;
		m_protein_range.end = m_protein_range.start + protein_size - 1;
	}
	
	all_size = rna_size + dna_size + protein_size;
	
}

void Structures::readPDB_rna(const std::filesystem::path& inputPath)
{
	ReadPDB inputPDB;
	inputPDB.read_file(inputPath);
	m_vAtoms_rna = inputPDB.m_vAtoms; 
	m_rna_seq = extract_rna_seq();
}

void Structures::readPDB_dna(const std::filesystem::path& inputPath)
{
	ReadPDB inputPDB;
	inputPDB.read_file(inputPath);
	m_vAtoms_dna = inputPDB.m_vAtoms;
	m_dna_seq = extract_dna_seq();
}

void Structures::readPDB_protein(const std::filesystem::path& inputPath)
{
	ReadPDB inputPDB;
	inputPDB.read_file(inputPath);
	m_vAtoms_protein = inputPDB.m_vAtoms;
	m_protein_seq = extract_protein_seq();
}

void Structures::sepRNA()
{
	m_vAtoms_rna.clear();
	m_vAtoms_rna.reserve(m_vAtoms_allStructures.size());
	for(auto const& atom : m_vAtoms_allStructures)
	{
		auto it = std::find(RNA.begin(), RNA.end(), del_spaces(atom.nuclName));
		if(it != RNA.end())
		{
			m_vAtoms_rna.emplace_back(atom);
		}
	}
}

void Structures::sepDNA()
{
	m_vAtoms_dna.clear();
	m_vAtoms_dna.reserve(m_vAtoms_allStructures.size());
	for(auto const& atom : m_vAtoms_allStructures)
	{
		auto it = std::find(DNA.begin(), DNA.end(), del_spaces(atom.nuclName));
		if(it != DNA.end())
		{
			m_vAtoms_dna.emplace_back(atom);
		}
	}
}

void Structures::sepProtein()
{
	m_vAtoms_protein.clear();
	m_vAtoms_protein.reserve(m_vAtoms_allStructures.size());
	for(auto const& atom : m_vAtoms_allStructures)
	{
		auto it = std::find(PROTEIN.begin(), PROTEIN.end(), del_spaces(atom.nuclName));
		if(it != PROTEIN.end())
		{
			m_vAtoms_protein.emplace_back(atom);
		}
	}
}

//void Structures::sepLigand()
//{ // it is not used, should be checked for the next versions
	//m_vAtoms_ligand.reserve(m_vAtoms_allStructures.size());
	//for(auto const& atom : m_vAtoms_allStructures)
	//{
		//auto itR = std::find(RNA.begin(), RNA.end(), del_spaces(atom.nuclName));
		//auto itD = std::find(DNA.begin(), DNA.end(), del_spaces(atom.nuclName));
		//auto itP = std::find(PROTEIN.begin(), PROTEIN.end(), del_spaces(atom.nuclName));
		//if(itR != RNA.end() && itD != DNA.end() && itP != PROTEIN.end())
		//{
			//m_vAtoms_ligand.emplace_back(atom);
		//}
	//}
	
	//for(const auto& atom : m_vAtoms_ligand)
	//{
		//cout << atom.chainId << " " << atom.atomNumber << " " << atom.nuclName << " " << atom.nuclNumber << endl;
	//}
	
//}

std::vector<std::string> Structures::extract_rna_seq()
{
	std::vector<std::string> rna_seq;
	std::string chain_seq;
	
	if(m_vAtoms_rna.size() == 0)
	{
		return rna_seq;
	}
	
	for(size_t i_atom  = 0; i_atom < m_vAtoms_rna.size()-1; ++i_atom)
	{
		if(m_vAtoms_rna[i_atom].nuclNumber != m_vAtoms_rna[i_atom + 1].nuclNumber)
		{
			chain_seq += del_spaces(m_vAtoms_rna[i_atom].nuclName);
			if(m_vAtoms_rna[i_atom].chainId != m_vAtoms_rna[i_atom + 1].chainId)
			{
				rna_seq.emplace_back(chain_seq);
				rna_size += chain_seq.size();
				chain_seq.clear();
			}
		}
	}
	chain_seq += del_spaces(m_vAtoms_rna[m_vAtoms_rna.size()-1].nuclName);
	rna_size += chain_seq.size();
	rna_seq.emplace_back(chain_seq);
	return rna_seq;
}

std::vector<std::string> Structures::extract_dna_seq()
{
	std::vector<std::string> dna_seq;
	std::string chain_seq;
	
	if(m_vAtoms_dna.size() == 0)
	{
		return dna_seq;
	}
	
	for(size_t i_atom  = 0; i_atom < m_vAtoms_dna.size()-1; ++i_atom)
	{
		if(m_vAtoms_dna[i_atom].nuclNumber != m_vAtoms_dna[i_atom + 1].nuclNumber)
		{
			chain_seq += DNA_SEQ[del_spaces(m_vAtoms_dna[i_atom].nuclName)];
			if(m_vAtoms_dna[i_atom].chainId != m_vAtoms_dna[i_atom + 1].chainId)
			{
				dna_seq.emplace_back(chain_seq);
				dna_size += chain_seq.size();
				chain_seq.clear();
			}
		}
	}
	chain_seq += DNA_SEQ[del_spaces(m_vAtoms_dna[m_vAtoms_dna.size()-1].nuclName)];
	dna_size += chain_seq.size();
	dna_seq.emplace_back(chain_seq);
	
	return dna_seq;
}

std::vector<std::string> Structures::extract_protein_seq()
{
	std::vector<std::string> protein_seq;
	std::string chain_seq;
	
	if(m_vAtoms_protein.size() == 0)
	{
		return protein_seq;
	}
	
	for(size_t i_atom  = 0; i_atom < m_vAtoms_protein.size()-1; ++i_atom)
	{
		if(m_vAtoms_protein[i_atom].nuclNumber != m_vAtoms_protein[i_atom + 1].nuclNumber)
		{
			chain_seq += PROTEIN_SEQ[del_spaces(m_vAtoms_protein[i_atom].nuclName)];
			if(m_vAtoms_protein[i_atom].chainId != m_vAtoms_protein[i_atom + 1].chainId)
			{
				protein_seq.emplace_back(chain_seq);
				protein_size += chain_seq.size();
				chain_seq.clear();
			}
		}
	}
	chain_seq += PROTEIN_SEQ[del_spaces(m_vAtoms_protein[m_vAtoms_protein.size()-1].nuclName)];
	protein_size += chain_seq.size();
	protein_seq.emplace_back(chain_seq);
	
	return protein_seq;
}

void Structures::read_restraint_file(const std::filesystem::path& restrainPath)
{
	std::ifstream res_file(restrainPath.string());
	std::string line;
	std::getline(res_file, line);
	
	//separate chains
	std::vector<std::string> v_chains_restraint;
	std::istringstream iss(line);
	std::string chain_restraint;
	while(iss >> chain_restraint)
	{
		v_chains_restraint.emplace_back(chain_restraint);
	}
	
	//// concatenate all sequence
	//std::string seq;
	//std::ostringstream oss;
	//std::copy(m_sequence_sep_chains.begin(), m_sequence_sep_chains.end()-1, std::ostream_iterator<std::string>(oss, " "));
	//oss << m_sequence_sep_chains[m_sequence_sep_chains.size()-1];
	//seq = oss.str();
	
	// check the constency of the restraint file and pdb file
	if(v_chains_restraint.size() != m_sequence_sep_chains.size())
	{
		std::cout << "restraint file and pdb file are not consistant\n";
		std::cout << "sequence : ";
		std::copy(m_sequence_sep_chains.begin(), m_sequence_sep_chains.end(), std::ostream_iterator<std::string>(std::cout, " "));
		std::cout << std::endl;
		std::cout << "restraint: ";
		std::copy(v_chains_restraint.begin(), v_chains_restraint.end(), std::ostream_iterator<std::string>(std::cout, " "));
		std::cout << std::endl;
		exit(EXIT_FAILURE);
	}
	
	// update the all strcutures vector
	std::vector<Atom> v_atoms_update;
	for(int i {0}; i < std::ssize(m_vAtoms_allStructures); ++i)
	{
		if(line[m_vAtoms_allStructures[i].nuclNumber - 1] == 'x') v_atoms_update.emplace_back( m_vAtoms_allStructures[i]);
	}
	
	m_vAtoms_allStructures = v_atoms_update;

		std::vector<Atom> v_atoms;
	std::string chain_seq;
	Residue residue;
	Chain chain;
	for(size_t i_atom  = 0; i_atom < m_vAtoms_allStructures.size()-1; ++i_atom)
	{
		v_atoms.emplace_back(m_vAtoms_allStructures[i_atom]);
		if(m_vAtoms_allStructures[i_atom].nuclNumber != m_vAtoms_allStructures[i_atom + 1].nuclNumber)
		{
			residue.v_atoms_in_residue.insert(residue.v_atoms_in_residue.end(), v_atoms.begin(), v_atoms.end());
			v_atoms.clear();
			residue.residue_name = del_spaces(m_vAtoms_allStructures[i_atom].nuclName);
			chain.v_residues_in_chain.emplace_back(residue);
			residue.v_atoms_in_residue.clear();
			if(m_vAtoms_allStructures[i_atom + 1].nuclNumber - m_vAtoms_allStructures[i_atom].nuclNumber > 1)
			{
				m_vChains.push_back(chain);
				chain.v_residues_in_chain.clear();
			}
		}
	}
	
	residue.v_atoms_in_residue.insert(residue.v_atoms_in_residue.end(), v_atoms.begin(), v_atoms.end());
	chain.v_residues_in_chain.emplace_back(residue);
	m_vChains.push_back(chain);

	std::cout << m_vChains.size() << std::endl;
	
	//update sequence
	std::vector<std::string> v_chains_seq_update;
	v_chains_seq_update.resize(m_sequence_sep_chains.size());
	for(int i {0}; i < std::ssize(v_chains_restraint); ++i)
	{
		for(int j {0}; j < std::ssize(v_chains_restraint[i]); ++j)
		{
			if(v_chains_restraint[i][j] == 'x') v_chains_seq_update[i] += m_sequence_sep_chains[i][j];
		}
	}
	m_sequence_sep_chains = v_chains_seq_update;
	//std::cout << "Udated sequence\n";
	//std::copy(m_sequence_sep_chains.begin(), m_sequence_sep_chains.end(), std::ostream_iterator<std::string>(std::cout, " "));
	//std::cout << '\n';
}

void Structures::read_contacts_interst(std::filesystem:: path const& ci_path)
{
	std::ifstream ci_file(ci_path.string());
	std::string line;
	
	while(std::getline(ci_file, line))
	{
		std::string s_num;
		std::istringstream iss{line};
		std::vector<int> v_tmp;
		while(iss >> s_num)
		{
			v_tmp.push_back(std::stoi(s_num));
		}
		// std::cout << v_tmp[0] << "--" << v_tmp[1] << std::endl;
		m_vp_contacts_interest.push_back({v_tmp[0], v_tmp[1]});
	}

	ci_file.close();
}

void Structures::check_contacts_interest(std::filesystem::path const& ci_path, std::filesystem::path const& out_path)
{
	read_contacts_interst(ci_path);
	int number_atoms_in_contact = 0;
	std::vector<int> v_ints(m_vp_contacts_interest.size());
	int m = 0, n_ok = 0;
	int seq_size = m_vChains[0].v_residues_in_chain.size();
	std::vector<bool> v_bools(seq_size,false);
	
	for(auto const& pair : m_vp_contacts_interest)
	{
		//std::cout << "size: " <<  std::ssize(m_vChains[0].v_residues_in_chain[pair.first-1].v_atoms_in_residue) << std::endl;
		for(int i_atom {0}; i_atom < std::ssize(m_vChains[0].v_residues_in_chain[pair.first-1].v_atoms_in_residue); ++i_atom)
		{
			int n = 0;
			
			for(int j_atom {0}; j_atom < std::ssize(m_vChains[0].v_residues_in_chain[pair.second-1].v_atoms_in_residue); ++j_atom)
			{
				number_atoms_in_contact = 0;
				auto dist = m_vChains[0].v_residues_in_chain[pair.first-1].v_atoms_in_residue[i_atom].coord.distance(
							m_vChains[0].v_residues_in_chain[pair.second-1].v_atoms_in_residue[j_atom].coord);
				if(dist < m_distance_threshold)
				{
					//std::cout << m_vChains[0].v_residues_in_chain[pair.first-1].residue_name << "--"
							  //<<pair.first  << "--" << i_atom << " " 
							  //<< m_vChains[0].v_residues_in_chain[pair.second-1].residue_name << "--"
							  //<< pair.second << "--" << j_atom << " --> " << dist << std::endl;
					++n;
					if(n == m_number_atoms_in_contact)
					{
						number_atoms_in_contact = n;
						break;
					}
				}
			}
			if(number_atoms_in_contact == m_number_atoms_in_contact)
			{
				// cout << pair.first << "<-->" << pair.second << " " << " " << distance << endl;
				v_ints[m] = 1;
				v_bools[pair.first-1] = true;
				v_bools[pair.second-1] = true;
				++n_ok;
				break;
			}
		}
		
		//auto dist =  m_vChains[0].v_residues_in_chain[pair.first-1].v_atoms_in_residue[pair.first].coord.distance(
							//m_vChains[0].v_residues_in_chain[pair.second-1].v_atoms_in_residue[pair.second].coord);
		//std::coou <
		//if(dist < m_distance_threshold)
		//{
			//std::cout << dist << " --> " << pair.first << "--" << pair.second << std::endl;
		//}
		++m;
	}
	
	std::ofstream out_file(out_path.string());
	std::ofstream out_file2(out_path.string() + "_2");
	
	std::copy(v_ints.begin(), v_ints.end() - 1, std::ostream_iterator<int>(out_file, ","));
	out_file << v_ints[v_ints.size() - 1] << std::endl;
	//out_file2
	std::copy(v_bools.begin(), v_bools.end() - 1, std::ostream_iterator<int>(out_file2, ","));
	out_file2 << v_bools[seq_size - 1] << std::endl;
	//out_file << n_ok << ",";
	//if(n_ok >= 4)
	//{
		//out_file << 1 << std::endl;
	//}
	//else
	//{
		//out_file << 0 << std::endl;
	//}

}

std::string Structures::get_requested_molecules() const
{
	return m_requested_molecules;
}

void Structures::set_requested_molecules(std::string requested_molecules)
{
	m_requested_molecules = requested_molecules;
}

double Structures::get_distance_threshold() const
{
	return m_distance_threshold;
}

void Structures::set_distance_threshold(double distance_threshold)
{
	m_distance_threshold = distance_threshold;
}


int Structures::get_number_atoms_in_contact() const
{
	return m_number_atoms_in_contact;
}

void Structures::set_number_atoms_in_contact(int number_atoms_in_contact)
{
	m_number_atoms_in_contact = number_atoms_in_contact;
}


size_t Structures::get_max_n_threads() const
{
	return m_max_n_threads;
}

void Structures::set_max_n_threads(size_t max_n_threads)
{
	m_max_n_threads = max_n_threads;
}

std::string Structures::get_calculation_type() const
{
	return m_calculation_type;
}

void Structures::set_calculation_type(std::string calculation_type)
{
	m_calculation_type = calculation_type;
}


bool Structures::get_is_plot_requested() const
{
	return m_is_plot_requested;
}

void Structures::set_is_plot_requested(bool is_plot_requested)
{
	m_is_plot_requested = is_plot_requested;
}

void Structures::create_all_seq()
{
	auto itR = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'R' || c == 'r');});

	if(itR != m_requested_molecules.end())
	{
		if(m_vAtoms_rna.size() > 0)
		{
			m_sequence_sep_chains.insert(m_sequence_sep_chains.end(), m_rna_seq.begin(), m_rna_seq.end());
		}
		else
		{
			cout << "Warning: There are no RNA in the input file\n";
		}
	}
	
	auto itD = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'D' || c == 'd');});
	if(itD != m_requested_molecules.end())
	{
		if(m_vAtoms_dna.size() > 0)
		{
			m_sequence_sep_chains.insert(m_sequence_sep_chains.end(), m_dna_seq.begin(), m_dna_seq.end());
		}
		else
		{
			cout << "Warning: There are no DNA in the input file\n";
		}
	}
	
	auto itP = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'P' || c == 'p');});
	if(itP != m_requested_molecules.end())
	{
		if(m_vAtoms_protein.size() > 0)
		{
			m_sequence_sep_chains.insert(m_sequence_sep_chains.end(), m_protein_seq.begin(), m_protein_seq.end());
		}
		else
		{
			cout << "Warning: There are no Protein in the input file\n";
		}
	}
}

void Structures::calc_distance_contact_all()
{
	initial_conatact_map(m_contact_map,all_size);

	size_t i,j; 
	i = 0;
	int number_atoms_in_contact = 0;
	for(size_t i_chain{ 0 }; i_chain < m_vChains.size(); ++i_chain)
	{
		for(size_t i_residue{ 0 }; i_residue < m_vChains[i_chain].v_residues_in_chain.size(); ++i_residue)
		{
			j = 0;
			for(size_t j_chain{ 0 }; j_chain  < m_vChains.size(); ++j_chain)
			{
				for(size_t j_residue{ 0 }; j_residue < m_vChains[j_chain].v_residues_in_chain.size(); ++j_residue)
				{
					number_atoms_in_contact = 0;
					for(size_t i_atom{ 0 }; i_atom < m_vChains[i_chain].v_residues_in_chain[i_residue].v_atoms_in_residue.size(); ++i_atom)
					{
						int n = 0;			
						for(size_t j_atom{ 0 }; j_atom < m_vChains[j_chain].v_residues_in_chain[j_residue].v_atoms_in_residue.size(); ++j_atom)
						{
							auto distance = m_vChains[i_chain].v_residues_in_chain[i_residue].v_atoms_in_residue[i_atom].coord.distance(
							m_vChains[j_chain].v_residues_in_chain[j_residue].v_atoms_in_residue[j_atom].coord);
							
					
							if(distance < m_distance_threshold)
							{
								++n;
								if(n == m_number_atoms_in_contact)
								{	
									number_atoms_in_contact = n; continue;
								}
							}
						}
					}
					if(number_atoms_in_contact == m_number_atoms_in_contact)
					{
						if(i == j)	m_contact_map[i][j] = m_contact_map[j][i] = 0;
						else if( abs(int(i - j)) == 1) m_contact_map[i][j] = m_contact_map[j][i] = 2; //immediate residues.
						else
						{
							m_contact_map[i][j] = m_contact_map[j][i] = 1;
							////rna-rna
							//if(i < m_rna_range.end && j < m_rna_range.end)
							//{
								//m_contact_map[i][j] = m_contact_map[j][i] = 1;
							//}
							
							////dna-dna
							//if(m_dna_range.start <= i && i < m_dna_range.end  && m_dna_range.start <= j && j < m_dna_range.end)
							//{
								//m_contact_map[i][j] = m_contact_map[j][i] = 1;
							//}
							
							////protein-protein
							//if(m_protein_range.start <= i && i < m_protein_range.end  && m_protein_range.start <= j && j < m_protein_range.end)
							//{
								//m_contact_map[i][j] = m_contact_map[j][i] = 3;
							//}
							
							////rna-dna
							//if((m_rna_range.end > i && m_dna_range.start <= j && j < m_dna_range.end)
							//|| (m_rna_range.end > j && m_dna_range.start <= i && i < m_dna_range.end))
							//{
								//m_contact_map[i][j] = m_contact_map[j][i] = 4;
							//}
							
							////rna-protein
							//if((m_rna_range.end > i && m_protein_range.start <= j && j < m_protein_range.end)
							//|| (m_rna_range.end > j && m_protein_range.start <= i && i < m_protein_range.end))
							//{
								//m_contact_map[i][j] = m_contact_map[j][i] = 5;
							//}
							
							////dna-protein
							//if((m_dna_range.end > i && m_protein_range.start <= j && j < m_protein_range.end)
							//|| (m_dna_range.end > j && m_protein_range.start <= i && i < m_protein_range.end))
							//{
								//m_contact_map[i][j] = m_contact_map[j][i] = 6;
							//}
							
						}
					}
					++j;	
				}
			}
			++i;
		}
	}
	
	sepMolecules_map();
}

void Structures::calc_distance_contact_inter()
{
	initial_conatact_map(m_contact_map,all_size);

	size_t i,j; 
	i = 0;

	int number_atoms_in_contact = 0;
	for(size_t i_chain{ 0 }; i_chain < m_vChains.size(); ++i_chain)
	{
		for(size_t i_residue{ 0 }; i_residue < m_vChains[i_chain].v_residues_in_chain.size(); ++i_residue)
		{
			j = 0;
			for(size_t j_chain{ 0 }; j_chain  < m_vChains.size(); ++j_chain)
			{
				for(size_t j_residue{ 0 }; j_residue < m_vChains[j_chain].v_residues_in_chain.size(); ++j_residue)
				{
					number_atoms_in_contact = 0;
					for(size_t i_atom{ 0 }; i_atom < m_vChains[i_chain].v_residues_in_chain[i_residue].v_atoms_in_residue.size(); ++i_atom)
					{
						int n = 0;			
						for(size_t j_atom{ 0 }; j_atom < m_vChains[j_chain].v_residues_in_chain[j_residue].v_atoms_in_residue.size(); ++j_atom)
						{							
							auto distance = m_vChains[i_chain].v_residues_in_chain[i_residue].v_atoms_in_residue[i_atom].coord.distance(
							m_vChains[j_chain].v_residues_in_chain[j_residue].v_atoms_in_residue[j_atom].coord);
							
							if(distance < m_distance_threshold)
							{
								++n;
								if(n == m_number_atoms_in_contact){ number_atoms_in_contact = n; continue;}
							}
						}
					}
					if(number_atoms_in_contact == m_number_atoms_in_contact)
					{
						if(i_chain == j_chain) m_contact_map[i][j] = m_contact_map[j][i] = 0;
						else
						{
							//rna-rna
							if(i < m_rna_range.end && j < m_rna_range.end)
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 1;
							}
							
							//dna-dna
							if(m_dna_range.start <= i && i < m_dna_range.end  && m_dna_range.start <= j && j < m_dna_range.end)
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 2;
							}
							
							//protein-protein
							if(m_protein_range.start <= i && i < m_protein_range.end  && m_protein_range.start <= j && j < m_protein_range.end)
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 3;
							}
							
							//rna-dna
							if((m_rna_range.end > i && m_dna_range.start <= j && j < m_dna_range.end)
							|| (m_rna_range.end > j && m_dna_range.start <= i && i < m_dna_range.end))
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 4;
							}
							
							//rna-protein
							if((m_rna_range.end > i && m_protein_range.start <= j && j < m_protein_range.end)
							|| (m_rna_range.end > j && m_protein_range.start <= i && i < m_protein_range.end))
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 5;
							}
							
							//dna-protein
							if((m_dna_range.end > i && m_protein_range.start <= j && j < m_protein_range.end)
							|| (m_dna_range.end > j && m_protein_range.start <= i && i < m_protein_range.end))
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 6;
							}
						}
					}
					++j;
				}
			}
			++i;
		}
	}

	sepMolecules_map();
}

void Structures::calc_distance_contact_intra()
{
	initial_conatact_map(m_contact_map,all_size);

	size_t i,j; 
	i = 0;

	int number_atoms_in_contact = 0;
	for(size_t i_chain{ 0 }; i_chain < m_vChains.size(); ++i_chain)
	{
		for(size_t i_residue{ 0 }; i_residue < m_vChains[i_chain].v_residues_in_chain.size(); ++i_residue)
		{
			j = 0;
			for(size_t j_chain{ 0 }; j_chain  < m_vChains.size(); ++j_chain)
			{
				for(size_t j_residue{ 0 }; j_residue < m_vChains[j_chain].v_residues_in_chain.size(); ++j_residue)
				{
					number_atoms_in_contact = 0;
					for(size_t i_atom{ 0 }; i_atom < m_vChains[i_chain].v_residues_in_chain[i_residue].v_atoms_in_residue.size(); ++i_atom)
					{
						if(i_chain != j_chain) continue;
						int n = 0;			
						for(size_t j_atom{ 0 }; j_atom < m_vChains[j_chain].v_residues_in_chain[j_residue].v_atoms_in_residue.size(); ++j_atom)
						{							
							auto distance = m_vChains[i_chain].v_residues_in_chain[i_residue].v_atoms_in_residue[i_atom].coord.distance(
							m_vChains[j_chain].v_residues_in_chain[j_residue].v_atoms_in_residue[j_atom].coord);
							
							if(distance < m_distance_threshold)
							{
								++n;
								if(n == m_number_atoms_in_contact){ number_atoms_in_contact = n; continue;}
							}
						}
					}
					if(number_atoms_in_contact == m_number_atoms_in_contact)
					{
						if(i == j)	m_contact_map[i][j]= m_contact_map[j][i] = 0;
						else if( abs(int(i - j)) == 1) m_contact_map[i][j] = m_contact_map[j][i] = 9; //immediate residues.
						else
						{
							//rna-rna
							if(i < m_rna_range.end && j < m_rna_range.end)
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 1;
							}
							
							//dna-dna
							if(m_dna_range.start <= i && i < m_dna_range.end  && m_dna_range.start <= j && j < m_dna_range.end)
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 2;
							}
							
							//protein-protein
							if(m_protein_range.start <= i && i < m_protein_range.end  && m_protein_range.start <= j && j < m_protein_range.end)
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 3;
							}
							
							//rna-dna
							if((m_rna_range.end > i && m_dna_range.start <= j && j < m_dna_range.end)
							|| (m_rna_range.end > j && m_dna_range.start <= i && i < m_dna_range.end))
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 4;
							}
							
							//rna-protein
							if((m_rna_range.end > i && m_protein_range.start <= j && j < m_protein_range.end)
							|| (m_rna_range.end > j && m_protein_range.start <= i && i < m_protein_range.end))
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 5;
							}
							
							//dna-protein
							if((m_dna_range.end > i && m_protein_range.start <= j && j < m_protein_range.end)
							|| (m_dna_range.end > j && m_protein_range.start <= i && i < m_protein_range.end))
							{
								m_contact_map[i][j] = m_contact_map[j][i] = 6;
							}
						}
						//m_contact_map[j][i] = 1;
					}
					++j;
				}
			}
			++i;
		}
	}
	
	sepMolecules_map();
}


void Structures::calc_distance_contact()
{
	if(m_calculation_type == "ALL") calc_distance_contact_all();
	if(m_calculation_type == "INTER") calc_distance_contact_inter();
	if(m_calculation_type == "INTRA") calc_distance_contact_intra();
}

void Structures::sepRNA_map()
{
	initial_conatact_map(m_contact_map_rna, rna_size);
	for(size_t i_res {m_rna_range.start}; i_res <= m_rna_range.end; ++i_res)
	{
		for(size_t j_res {i_res + 1}; j_res <= m_rna_range.end; ++j_res)
		{
			m_contact_map_rna[i_res][j_res] = m_contact_map[i_res][j_res];
			m_contact_map_rna[j_res][i_res] = m_contact_map[j_res][i_res];
		}
	}
}

void Structures::sepDNA_map()
{
	initial_conatact_map(m_contact_map_dna, dna_size);
	for(size_t i_res {m_dna_range.start}; i_res <= m_dna_range.end; ++i_res)
	{
		int i = i_res - m_rna_range.end;
		for(size_t j_res {i_res + 1}; j_res <= m_dna_range.end; ++j_res)
		{
			int j = j_res - m_rna_range.end;
			m_contact_map_dna[i][j] = m_contact_map[i_res][j_res];
			m_contact_map_dna[j][i] = m_contact_map[j_res][i_res];
		}
	}
}

void Structures::sepProtein_map()
{
	initial_conatact_map(m_contact_map_protein, protein_size);
	for(size_t i_res {m_protein_range.start}; i_res < m_protein_range.end; ++i_res)
	{
		int i = i_res - m_rna_range.end - m_dna_range.end;
		for(size_t j_res {i_res + 1}; j_res < m_protein_range.end; ++j_res)
		{
			int j = j_res - m_rna_range.end - m_dna_range.end;
			m_contact_map_protein[i][j] = m_contact_map[i_res][j_res];
			m_contact_map_protein[j][i] = m_contact_map[j_res][i_res];
		}
	}
}

void Structures::sepRNADNA_map()
{
	initial_conatact_map(m_contact_map_rna_dna, all_size);
	for(size_t i_res {0}; i_res < all_size; ++i_res)
	{
		for(size_t j_res {i_res + 1}; j_res < all_size; ++j_res)
		{
			if(i_res < m_rna_range.end && m_dna_range.start < j_res && j_res < m_protein_range.start)
			{
				m_contact_map_rna_dna[i_res][j_res] = m_contact_map[i_res][j_res];
				m_contact_map_rna_dna[j_res][i_res] = m_contact_map[j_res][i_res];
			}
			else
			{
				m_contact_map_rna_dna[i_res][j_res] = 0;
				m_contact_map_rna_dna[j_res][i_res] = 0;
			}
		}
	}
}

void Structures::sepRNAProtein_map()
{
	initial_conatact_map(m_contact_map_rna_protein, all_size);
	for(size_t i_res {0}; i_res < all_size; ++i_res)
	{
		for(size_t j_res {i_res + 1}; j_res < all_size; ++j_res)
		{
			if(i_res < m_rna_range.end && j_res >= m_protein_range.start)
			{
				m_contact_map_rna_protein[i_res][j_res] = m_contact_map[i_res][j_res];
				m_contact_map_rna_protein[j_res][i_res] = m_contact_map[j_res][i_res];
			}
			else
			{
				m_contact_map_rna_protein[i_res][j_res] = 0;
				m_contact_map_rna_protein[j_res][i_res] = 0;
			}
		}
	}
}

void Structures::sepDNAProtein_map()
{
	initial_conatact_map(m_contact_map_dna_protein, all_size);
	for(size_t i_res {0}; i_res < all_size; ++i_res)
	{
		for(size_t j_res {i_res + 1}; j_res < all_size; ++j_res)
		{
			if(m_rna_range.end < i_res && i_res < m_dna_range.end && j_res > m_protein_range.start)
			{
				m_contact_map_dna_protein[i_res][j_res] = m_contact_map[i_res][j_res];
				m_contact_map_dna_protein[j_res][i_res] = m_contact_map[j_res][i_res];
			}
			else
			{
				m_contact_map_dna_protein[i_res][j_res] = 0;
				m_contact_map_dna_protein[j_res][i_res] = 0;
			}
		}
	}
}


void Structures::sepMolecules_map()
{	
	sepRNA_map();
	sepRNADNA_map();
	sepRNAProtein_map();
	
	sepDNA_map();
	sepDNAProtein_map();
	
	sepProtein_map();
}

std::vector<std::string> Structures::make_raw_structure(const std::vector<std::string>& v_seq) const
{
	std::vector<std::string> v_raw_structure;
	for(auto const& chain_seq : v_seq)
	{
		std::string raw_structure(chain_seq.size(), '.');
		v_raw_structure.emplace_back(raw_structure);
	}
	return v_raw_structure;
}



std::vector<std::string> Structures::parse_mole_binary(const std::vector<std::vector<int>>& contact_map, const std::vector<std::string>& v_seq) const
{	
	std::string raw_struct(contact_map.size(), '.');
	for(size_t i { 0 }; i < contact_map.size(); ++i)
	{
		for(size_t j { i+1 }; j < contact_map.size(); ++j)
		{
			if(contact_map[i][j] > 0 && contact_map[i][j] < 9)
			{
				raw_struct[i] = 'X';
				raw_struct[j] = 'X';
			}
		}
	}
	std::vector<std::string> v_binary = structure_string2vec(raw_struct, v_seq);
	return v_binary;
}

void Structures::write_binary_format(const std::filesystem::path& outputPath) const
{
	std::vector<std::string> v_binary;
	auto itR = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'R' || c == 'r');});
	auto itD = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'D' || c == 'd');});
	auto itP = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'P' || c == 'p');});
	
	if(itR != m_requested_molecules.end())
	{
		auto v_rna_binary = parse_mole_binary(m_contact_map_rna, m_rna_seq);
		if(v_rna_binary.size() > 0)
		{
			std::ofstream rna_outFile(outputPath.string() + "_rna_" + m_calculation_type + ".xo");
			std::copy(v_rna_binary.begin(), v_rna_binary.end() - 1, std::ostream_iterator<std::string>(rna_outFile, " "));
			rna_outFile << v_rna_binary[v_rna_binary.size() - 1];
			rna_outFile.close();
		}
		else
		{
			cout << "Warning: your request is aborted\n";
			cout << "There is no RNA in the input pdb file\n";
		}
	}
	
	if(itD != m_requested_molecules.end())
	{
		auto v_dna_binary = parse_mole_binary(m_contact_map_dna, m_dna_seq);
		if(v_dna_binary.size() > 0)
		{
			std::ofstream dna_outFile(outputPath.string() + "_dna_" + m_calculation_type + ".xo");
			std::copy(v_dna_binary.begin(), v_dna_binary.end() - 1, std::ostream_iterator<std::string>(dna_outFile, " "));
			dna_outFile << v_dna_binary[v_dna_binary.size() - 1];
			dna_outFile.close();
		}
		else
		{
			cout << "Warning: your request is aborted\n";
			cout << "There is no DNA in the input pdb file\n";
		}
	}
	
	if(itP != m_requested_molecules.end())
	{
		auto v_protein_binary = parse_mole_binary(m_contact_map_protein, m_protein_seq);
		if(v_protein_binary.size() > 0)
		{
			std::ofstream protein_outFile(outputPath.string() + "_protein_" + m_calculation_type + ".xo");
			std::copy(v_protein_binary.begin(), v_protein_binary.end() - 1, std::ostream_iterator<std::string>(protein_outFile, " "));
			protein_outFile << v_protein_binary[v_protein_binary.size() - 1];
			protein_outFile.close();
		}
		else
		{
			cout << "Warning: your request is aborted\n";
			cout << "There is no protein in the input pdb file\n";
		}
	}
	
	
	
	if(itR != m_requested_molecules.end() && itD != m_requested_molecules.end())
	{
		auto v_rna_dna_binary = parse_mole_binary(m_contact_map_rna_dna, m_sequence_sep_chains);
		if(v_rna_dna_binary.size() > 0)
		{
			std::ofstream rna_dna_outFile(outputPath.string() + "_rna_dna_" + m_calculation_type + ".xo");
			std::copy(v_rna_dna_binary.begin(), v_rna_dna_binary.end() - 1, std::ostream_iterator<std::string>(rna_dna_outFile, " "));
			rna_dna_outFile << v_rna_dna_binary[v_rna_dna_binary.size() - 1];
			rna_dna_outFile.close();
		}
	}
	
	if(itR != m_requested_molecules.end() && itP != m_requested_molecules.end())
	{
		auto v_rna_protein_binary = parse_mole_binary(m_contact_map_rna_protein, m_sequence_sep_chains);
		if(v_rna_protein_binary.size() > 0)
		{
			std::ofstream rna_protein_outFile(outputPath.string() + "_rna_protein_" + m_calculation_type + ".xo");
			std::copy(v_rna_protein_binary.begin(), v_rna_protein_binary.end() - 1, std::ostream_iterator<std::string>(rna_protein_outFile, " "));
			rna_protein_outFile << v_rna_protein_binary[v_rna_protein_binary.size() - 1];
			rna_protein_outFile.close();
		}
	}
	
		
	if(itD != m_requested_molecules.end() && itP != m_requested_molecules.end())
	{
		auto v_dna_protein_binary = parse_mole_binary(m_contact_map_dna_protein, m_sequence_sep_chains);
		if(v_dna_protein_binary.size() > 0)
		{
			std::ofstream dna_protein_outFile(outputPath.string() + "_dna_protein_" + m_calculation_type + ".xo");
			std::copy(v_dna_protein_binary.begin(), v_dna_protein_binary.end() - 1, std::ostream_iterator<std::string>(dna_protein_outFile, " "));
			dna_protein_outFile << v_dna_protein_binary[v_dna_protein_binary.size() - 1];
			dna_protein_outFile.close();
		}
	}
	
	
	std::string seq_filename = outputPath.string() + ".seq";
	write_seq(seq_filename);
}


void Structures::write_as_dot_bracket(const std::filesystem::path& outputPath) const
{
	std::vector<std::string> v_dot_barckets;
	
	auto itR = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'R' || c == 'r');});
	auto itD = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'D' || c == 'd');});
	auto itP = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'P' || c == 'p');});
	
	if(itR != m_requested_molecules.end() && itP != m_requested_molecules.end() && itD == m_requested_molecules.end())
	{
		std::ofstream outFile(outputPath.string() + "_rna_protein.dbn");
		for(size_t i { 0 }; i < m_contact_map_rna_protein.size(); ++i)
		{
			for(size_t j { i+1 }; j < m_contact_map_rna_protein.size(); ++j)
			{
				//std::vector<std::string> v_dot_bracket_rna_protein;
				if(m_contact_map_rna_protein[i][j] > 1 && m_contact_map[i][j] < 9)
				{
					std::string raw_struct(m_contact_map_rna_protein.size(), '.');
					raw_struct[i] = '(';
					raw_struct[j] = ')';
					auto tmp_v = structure_string2vec(raw_struct, m_sequence_sep_chains);
					std::copy(tmp_v.begin(), tmp_v.end()-1, std::ostream_iterator<std::string>(outFile, " "));
					outFile << tmp_v[tmp_v.size() - 1] << "\n";
				}
			}
		}
		outFile.close();
	}
	
		
	if(itR == m_requested_molecules.end() && itP != m_requested_molecules.end() && itD != m_requested_molecules.end())
	{
		std::ofstream outFile(outputPath.string() + "_dna_protein.dbn");
		for(size_t i { 0 }; i < m_contact_map_dna_protein.size(); ++i)
		{
			for(size_t j { i+1 }; j < m_contact_map_dna_protein.size(); ++j)
			{
				if(m_contact_map_dna_protein[i][j] > 1 && m_contact_map[i][j] < 9)
				{
					std::string raw_struct(m_contact_map_dna_protein.size(), '.');
					raw_struct[i] = '(';
					raw_struct[j] = ')';
					auto tmp_v = structure_string2vec(raw_struct, m_sequence_sep_chains);
					std::copy(tmp_v.begin(), tmp_v.end()-1, std::ostream_iterator<std::string>(outFile, " "));
					outFile << tmp_v[tmp_v.size() - 1] << "\n";
				}
			}
		}
		outFile.close();
	}
	
	
	
		
	if(itR != m_requested_molecules.end() && itP == m_requested_molecules.end() && itD != m_requested_molecules.end())
	{
		std::ofstream outFile(outputPath.string() + "_rna_dna.dbn");
		for(size_t i { 0 }; i < m_contact_map_rna_dna.size(); ++i)
		{
			for(size_t j { i+1 }; j < m_contact_map_rna_dna.size(); ++j)
			{
				//std::vector<std::string> v_dot_bracket_rna_protein;
				if(m_contact_map_rna_dna[i][j] > 1 && m_contact_map[i][j] < 9)
				{
					std::string raw_struct(m_contact_map_rna_dna.size(), '.');
					raw_struct[i] = '(';
					raw_struct[j] = ')';
					auto tmp_v = structure_string2vec(raw_struct, m_sequence_sep_chains);
					std::copy(tmp_v.begin(), tmp_v.end()-1, std::ostream_iterator<std::string>(outFile, " "));
					outFile << tmp_v[tmp_v.size() - 1] << "\n";
				}
			}
		}
		outFile.close();
	}
	
	
	std::ofstream outFile(outputPath.string() + "_all.dbn");
	for(size_t i { 0 }; i < m_contact_map.size(); ++i)
	{
		for(size_t j { i+1 }; j < m_contact_map.size(); ++j)
		{
			//std::vector<std::string> v_dot_bracket_rna_protein;
			if(m_contact_map[i][j] > 1 && m_contact_map[i][j] < 9)
			{
				std::string raw_struct(m_contact_map.size(), '.');
				raw_struct[i] = '(';
				raw_struct[j] = ')';
				auto tmp_v = structure_string2vec(raw_struct, m_sequence_sep_chains);
					std::copy(tmp_v.begin(), tmp_v.end()-1, std::ostream_iterator<std::string>(outFile, " "));
					outFile << tmp_v[tmp_v.size() - 1] << "\n";
			}
		}
	}
	
	outFile.close();
}



void Structures::write_seq(const std::filesystem::path& outputPath) const
{
	std::ofstream outFile(outputPath.string());
	if(std::ssize(m_sequence_sep_chains) == 1)
	{
		outFile << m_sequence_sep_chains[m_sequence_sep_chains.size()-1]<<endl;
	}
	else
	{
		std::copy(m_sequence_sep_chains.begin(), m_sequence_sep_chains.end()-1, std::ostream_iterator<std::string>(outFile, " "));
		outFile << m_sequence_sep_chains[m_sequence_sep_chains.size()-1]<<endl;
	}
}


void Structures::write_map(const std::filesystem::path& outputPath, const std::vector<std::vector<int>>& contact_map) const
{
	std::ofstream outFile(outputPath.string());
	std::ofstream outFile2(outputPath.string() + "_2");
	for(size_t i { 0 }; i < contact_map.size(); ++i)
	{
		for(size_t j { 0 }; j < contact_map.size(); ++j)
		{
			outFile << contact_map[i][j] << " ";
			outFile2 << contact_map[i][j];
		}
		outFile << endl;
		outFile2 << endl;
	}
	
	outFile.close();
	outFile2.close();
}

void Structures::write_map(const std::filesystem::path& outputPath) const
{
	//auto itR = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'R' || c == 'r');});
	//auto itD = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'D' || c == 'd');});
	//auto itP = std::find_if(m_requested_molecules.begin(), m_requested_molecules.end(), [](char c){return (c == 'P' || c == 'p');});
	
	//if(itR != m_requested_molecules.end())
	//{
			//std::ostringstream oss_rna;
			//std::ostringstream oss_plot_rna;
			//oss_rna << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_RNA_.map";
			//oss_plot_rna << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_RNA_.png";
			//write_map(oss_rna.str(), m_contact_map_rna);
			//cout << "\n\t\t **** plot -R- contacts ****\n";
			//plot(oss_plot_rna.str(), oss_rna.str());
	//}
	
	//if(itD != m_requested_molecules.end())
	//{
		//std::ostringstream oss_dna;
		//std::ostringstream oss_plot_dna;
		//oss_dna << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_DNA.map";
		//oss_plot_dna << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_DNA.png";
		//write_map(oss_dna.str(), m_contact_map_dna);
		//cout << "\n\t\t **** plot -D- contacts ****\n";
		//plot(oss_plot_dna.str(), oss_dna.str());
	//}
	
	//if(itP != m_requested_molecules.end())
	//{
		//std::ostringstream oss_protein;
		//std::ostringstream oss_plot_protein;
		//oss_protein << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_Protein.map";
		//oss_plot_protein << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_Protein.png";
		//write_map(oss_protein.str(), m_contact_map_protein);
		//cout << "\n\t\t **** plot -P- contacts ****\n";
		//plot(oss_plot_protein.str(), oss_protein.str());
	//}
	
	//if(itR != m_requested_molecules.end() && itD != m_requested_molecules.end())
	//{
		//std::ostringstream oss_rna_dna;
		//std::ostringstream oss_plot_rna_dna;
		//oss_rna_dna << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_RNA_DNA.map";
		//oss_plot_rna_dna << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_RNA_DNA.png";
		//write_map(oss_rna_dna.str(), m_contact_map_rna_dna);
		//if(m_calculation_type == "INTER" || m_calculation_type == "ALL")
		//{
			//cout << "\n\t\t **** plot R-D contacts ****\n";
			//plot(oss_plot_rna_dna.str(), oss_rna_dna.str());
		//}
	//}
	
	//if(itR != m_requested_molecules.end() && itP != m_requested_molecules.end())
	//{
		//std::ostringstream oss_rna_protein;
		//std::ostringstream oss_plot_rna_protein;
		//oss_rna_protein << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_RNA_Protein.map";
		//oss_plot_rna_protein << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_RNA_Protein.png";
		//write_map(oss_rna_protein.str(), m_contact_map_rna_protein);
		//if(m_calculation_type == "INTER" || m_calculation_type == "ALL")
		//{
			//cout << "\n\t\t **** plot R-P contacts ****\n";
			//plot(oss_plot_rna_protein.str(), oss_rna_protein.str());
		//}
	//}
	
	//if(itD != m_requested_molecules.end() && itP != m_requested_molecules.end())
	//{
		//std::ostringstream oss_dna_protein;
		//std::ostringstream oss_plot_dna_protein;
		//oss_dna_protein << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_DNA_Protein.map";
		//oss_plot_dna_protein << outputPath.parent_path().string() << "/" << outputPath.stem().string() << "_DNA_Protein.png";
		//write_map(oss_dna_protein.str(), m_contact_map_dna_protein);
		//if(m_calculation_type == "INTER" || m_calculation_type == "ALL")
		//{
			//cout << "\n\t\t **** plot D-P contacts ****\n";
			//plot(oss_plot_dna_protein.str(), oss_dna_protein.str());
		//}
	//}
	
	std::ostringstream oss_all;
	std::ostringstream oss_plot_all;
	oss_all << outputPath.string() << "_All.map";
	oss_plot_all << outputPath.string() << "_All.png";
	write_map(oss_all.str(), m_contact_map);
	cout << "\n\t\t **** plot all contacts ****\n";
	plot(oss_plot_all.str(), oss_all.str());
	
	//auto v_binary = parse_mole_binary(m_contact_map_rna_protein, m_sequence_sep_chains);
}


void Structures::plot(const std::filesystem::path& outputPath, const std::filesystem::path& contact_map_path) const
{
	gnuplot gnu;
	std::string input = "splot \"" + contact_map_path.string() + "\" matrix with image pixels";
	std::string output = "set output \"" + outputPath.string() + "\"";
	
	
	gnu("set lmargin 0");
	gnu("set rmargin 0");
	gnu("set tmargin 0");
	gnu("set bmargin 0");
	gnu("set term png size 1000,1000");
	gnu(output);
	gnu("set pm3d map");
	gnu("set palette defined( 0 \"white\", 1 \"black\", 2 \"grey\")");
	gnu("set yrange [0:*] reverse");
	gnu("unset ytics");
	gnu("unset xtics");
	gnu("unset colorbox");
	//gnu("unset border");
	gnu(input);	
}
