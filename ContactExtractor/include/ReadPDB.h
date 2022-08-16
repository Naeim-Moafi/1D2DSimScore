 
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

s




#pragma once
#define MAX_N_ATOMS_IN_PDB 16000 //it is not limiting constant any more, since chainging static tables to C++ vectors

#include <filesystem>
#include <fstream>
#include <iostream>
#include "Point3D.h"
#include <string>
#include <vector>

using std::cout;
using std::endl;


//typedef struct
//{
  //int start_index;
  //int n_atoms;
//} Residue;

//typedef struct
//{
  //int start_index;
  //int n_residues,n_residues_sum;
//} Chain;


////typedef struct{
	////double x, y, z;
////} Point3D;

struct Range {
	size_t start;
	size_t end;
	Range() {start = end = 0;};
};

struct Atom{
//class represents one atom of any type
// class will surely grow in the future
public:
	Point3D coord, coordBackup;
	std::string atomName, nuclName, chainId;
	int atomNumber, nuclNumber;
	double occupancy, bfactor;
	
	bool operator==(const Atom& atom)
	{
		bool coord_bool = (this->coord.x == atom.coord.x && this->coord.y == atom.coord.y && this->coord.z == atom.coord.z);
		bool chainId_bool = (this->chainId == atom.chainId);
		bool names_bool = (this->atomName == atom.atomName && this->nuclName == atom.nuclName);
		bool num_bool = (this->atomNumber == atom.atomNumber && this->nuclNumber == atom.nuclNumber);
		return (chainId_bool && coord_bool && names_bool && num_bool);
	}
};



class ReadPDB
{
	private:
		int n_atoms,n_residues,n_chains,selected_chain;
		bool fill_atom_strc(Atom& atom, const std::string& line);
		//int split_into_names(char **v_names, const char *i_names);
	public:
		std::vector<Atom> m_vAtoms;
		std::vector<std::vector<Atom>> m_vAtomsInChains;
		void read_file(const std::filesystem::path& pdbPath);
		int get_n_atoms() const;
		//void display_residue_lines(int i_res);
		
		
		//char get_chain_id(int i_chain);
		
		//void test_split_string(char *inp_string);
		ReadPDB()=default;
		~ReadPDB()=default;
};

bool are_there_spaces(const std::string& inp_str);
std::string del_spaces(const std::string& inp_str);
int split_string(std::vector<std::string>& output_vector, const std::string& string2parse, char delim);
