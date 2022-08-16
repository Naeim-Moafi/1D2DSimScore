 
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
#define MAX_N_ATOMS_IN_PDB 16000 //it is not limiting constant any more, since chainging static tables to C++ vectors


#include <filesystem>
#include <fstream>
#include <iostream>
//#include "Point3D.h"
#include <string>
#include <vector>

using std::cout;
using std::endl;


typedef struct
{
  int start_index;
  int n_atoms;
} Residue;

typedef struct
{
  int start_index;
  int n_residues,n_residues_sum;
} Chain;


typedef struct{
	double x, y, z;
} Point3D;

struct Range {
	int start;
	int term;
	int chainNum;
	
	int startBase;
	int termBase;
};

typedef struct{
//class represents one atom of any type
// class will surely grow in the future
public:
	Point3D coord, coordBackup;
	std::string atomName, nuclName, chainId;
	int atomNumber, nuclNumber;
	double occupancy, bfactor;
}Atom;



class ReadPDB
{
	private:
		int n_atoms,n_residues,n_chains,selected_chain;
		bool fill_atom_strc(Atom& atom, const std::string& line);
		std::vector<Residue> m_vResidues;
		std::vector<Chain> m_vChains;
		//int split_into_names(char **v_names, const char *i_names);
	public:
		std::vector<Atom> m_vAtoms;
		std::vector<std::vector<Atom>> m_vAtomsInChains;
		void read_file(const std::filesystem::path& pdbPath);
		int get_n_atoms() const;
		int get_n_chains()const;
		int get_n_residues()const;
		int get_n_residues_in_chain(int i_chain) const;
		int select_chain_for_reading(int i_chain);
		int get_atom(Atom& i_atom, int i_res, const std::string& i_name)const;
		int get_atom(Atom& i_atom, int i_res, const std::string& i_name, int i_chain) const;
		int get_atom_from_selected_chain(Atom& i_atom, int i_res, const std::string& i_name) const;
		//void display_line(int i_line);
		//void display_residue_lines(int i_res);
		
		
		//char get_chain_id(int i_chain);
		
		//void test_split_string(char *inp_string);
		ReadPDB()=default;
		~ReadPDB()=default;
};

bool are_there_spaces(const std::string& inp_str);
std::string del_spaces(const std::string& inp_str);
int split_string(std::vector<std::string>& output_vector, const std::string& string2parse, char delim);
