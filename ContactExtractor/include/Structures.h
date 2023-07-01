/*
 * Author: S. Naeim Moafinejad
 */
 
  
/*
 * International Institute of Molecular and Cell Biology (IIMCB)
 * Copyright [2022] [IIMCB]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *        http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


#pragma once


#include "gnuplot.h"
#include <map>
#include "ReadPDB.h"
#include <tuple>
#include <utility>

// forward declaration for testing
class TestStructres;

struct Residue
{
	std::vector<Atom> v_atoms_in_residue;
	std::string residue_name;
};

struct Chain
{
	std::vector<Residue> v_residues_in_chain;
	std::string ID;
};

struct ContactMap
{
	int first_chain_no, second_chain_no, first_resNumber, second_resNumber;
};

class Structures
{
	public:
		bool update = false;
		std::string selection_path;
		std::vector<std::string> m_sequence_sep_chains;
		std::vector<std::string> m_rna_seq;
		std::vector<std::string> m_dna_seq;
		std::vector<std::string> m_protein_seq;
		std::vector<Atom> m_vAtoms_allStructures;
		std::vector<Atom> m_vAtoms_rna;
		std::vector<Atom> m_vAtoms_rna_update;
		std::vector<Atom> m_vAtoms_dna;
		std::vector<Atom> m_vAtoms_protein;
		std::vector<Atom> m_vAtoms_ligand;
		std::vector<std::vector<Atom>> m_tRNA_only;
		std::vector<std::vector<int>> m_contact_map;
		std::vector<std::vector<int>> m_contact_map_rna;
		std::vector<std::vector<int>> m_contact_map_dna;
		std::vector<std::vector<int>> m_contact_map_protein;
		std::vector<std::vector<int>> m_contact_map_rna_dna;
		std::vector<std::vector<int>> m_contact_map_rna_protein;
		std::vector<std::vector<int>> m_contact_map_dna_protein;
		
		
		[[nodiscard]]std::string get_requested_molecules() const;
		void set_requested_molecules(std::string requested_molecules);
		[[nodiscard]]double get_distance_threshold() const;
		void set_distance_threshold(double distance_threshold);
		[[nodiscard]]int get_number_atoms_in_contact() const;
		void set_number_atoms_in_contact(int number_atoms_in_contact);
		[[nodiscard]]std::string get_calculation_type() const;
		void set_calculation_type(std::string calculation_type);
		[[nodiscard]]bool get_is_plot_requested() const;
		void set_is_plot_requested(bool is_plot_requested);
		size_t get_max_n_threads() const;
		void set_max_n_threads(size_t max_n_threads);
		
		void readPDB_allStructures(const std::filesystem::path& inputPath);
		void readPDB_rna(const std::filesystem::path& inputPath);
		void readPDB_dna(const std::filesystem::path& inputPath);
		void readPDB_protein(const std::filesystem::path& inputPath);
		void readPDB_ligand(const std::filesystem::path& inputPath);
		
		void read_restraint_file(const std::filesystem::path& inputPath);
		void check_contacts_interest(std::filesystem::path const& ci_path, std::filesystem::path const& out_path);
		
		void calc_distance_contact_all();
		void calc_distance_contact_inter();
		void calc_distance_contact_intra();
		void calc_distance_contact();
		[[nodiscard]]std::vector<ContactMap> 
		calc_distance_contact(const std::vector<Chain>& v_chains_first, const std::vector<Chain>& v_chains_second) const;
		void write_map(const std::filesystem::path& ouputpath, const std::vector<std::vector<int>>& contact_map) const;
		void write_map(const std::filesystem::path& ouputpath) const;
		void write_binary_format(const std::filesystem::path& outputPath) const;
		void write_as_dot_bracket(const std::filesystem::path& outputPath) const;
		void write_seq(const std::filesystem::path& outputPath) const;
		void create_all_seq();
		
		Structures() = default;
		~Structures() = default;
	#ifdef __TEST__
	public:
	else:
	private:
	#endif
		bool m_is_plot_requested = false;
		std::string m_requested_molecules;
		std::string m_calculation_type;
		double m_distance_threshold;
		int m_number_atoms_in_contact;
		size_t rna_size=0, dna_size=0, protein_size=0, all_size=0;
		Range m_rna_range, m_dna_range, m_protein_range;
		size_t m_max_n_threads;
		std::vector<Chain> m_vChains_rna;
		std::vector<Chain> m_vChains_dna;
		std::vector<Chain> m_vChains_protein;
		std::vector<Chain> m_vChains;
		std::vector<pair<int, int>> m_vp_contacts_interest;
	
		void sepRNA();
		void sepDNA();
		void sepProtein();
		void sepLigand();
		void sepRNA_map();
		void sepDNA_map();
		void sepProtein_map();
		void sepRNADNA_map();
		void sepRNAProtein_map();
		void sepDNAProtein_map();
		void sepMolecules_map();
		void read_contacts_interst(std::filesystem::path const& ci_path);
		void plot(const std::filesystem::path& outputPath, const std::filesystem::path& contact_map_path) const;
		std::vector<std::string> extract_rna_seq();
		std::vector<std::string> extract_dna_seq();
		std::vector<std::string> extract_protein_seq();
		std::vector<std::string> extract_ligand_seq();
		std::vector<std::string> parse_mole_binary(const std::vector<std::vector<int>>& contact_map, const std::vector<std::string>& v_seq) const;
		std::vector<std::string> make_raw_structure(const std::vector<std::string>& v_seq) const;
		friend class TestStructres;
		
};

void initial_conatact_map(std::vector<std::vector<int>>& contact_map, size_t size);
void initial_conatact_map(std::vector<std::vector<int>>& contact_map, size_t size1, size_t size2);
std::vector<std::string> structure_string2vec(const std::string& inp_struct, const std::vector<std::string>& v_seq);
void clear_v_string(std::vector<std::string>& v_str);