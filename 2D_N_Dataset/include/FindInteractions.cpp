#include "FindInteractions.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <cstdio>

using std::cout;
using std::endl;
using std::cerr;


const std::string closeSymbols = ")]}>abcdefghijklmnopqrstuvwxyz";
const std::string openSymbols = "([{<ABCDEFGHIGKLMNOPQRSTUVWXYZ";
const std::string possibleChainID = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
const std::vector<std::string> vBPh = {"W_0BPh", "W_1BPh", "W_2BPh", "W_345BPh", "W_6BPh", "W_789BPh",
									   "S_0BPh", "S_1BPh", "S_2BPh", "S_345BPh", "S_6BPh", "S_789BPh",
									   "H_0BPh", "H_1BPh", "H_2BPh", "H_345BPh", "H_6BPh", "H_789BPh"};
const std::vector<std::string> vBR = {"W_0BR", "W_1BR", "W_2BR", "W_345BR", "W_6BR", "W_789BR",
									   "S_0BR", "S_1BR", "S_2BR", "S_345BR", "S_6BR", "S_789BR",
									   "H_0BR", "H_1BR", "H_2BR", "H_345BR", "H_6BR", "H_789BR"};
									   
									   
														 

//=======================================================================================================
CLARNA_all::FindInteractions::FindInteractions()
	: CLARNA::FindInteractions::FindInteractions(){}

CLARNA_all::FindInteractions::FindInteractions(bool isWobble_canonical)
	: CLARNA::FindInteractions::FindInteractions(isWobble_canonical){}


bool CLARNA_all::FindInteractions::get_isWobble_canonical() const
{
	return CLARNA::FindInteractions::get_isWobble_canonical();
}

void CLARNA_all::FindInteractions::set_isWobble_canonical(bool isWobble_canonical)
{
	CLARNA::FindInteractions::set_isWobble_canonical(isWobble_canonical);
}

int CLARNA_all::FindInteractions::get_number_involved_faces_edges() const
{
	return CLARNA::FindInteractions::get_number_involved_faces_edges();
}

void CLARNA_all::FindInteractions::set_number_involved_faces_edges(int number_involved_faces_edges)
{
	CLARNA::FindInteractions::set_number_involved_faces_edges(number_involved_faces_edges);
}

std::string CLARNA_all::FindInteractions::get_extension() const
{
	return m_extension;
}

void CLARNA_all::FindInteractions::set_extension(const std::string& extension)
{
	m_extension = extension;
}


//void CLARNA_all::FindInteractions::init_sequence(const std::filesystem::path& path)
//{
	//if(std::filesystem::is_directory(path))
	//{
		////cout << "finding sequence file in folder and read it\n";
		//const std::filesystem::path target_dir { path };
		//// iterate all over the target directory to find the sequence file
		//for(const auto& dir_entry : std::filesystem::directory_iterator { target_dir })
		//{
			//const std::filesystem::path tmp_path { dir_entry };
			//if(tmp_path.extension() == ".pdb")
			//{
				//[[maybe_unused]] auto returnValue = CLARNA::FindInteractions::readPDBFile(tmp_path);
				//m_sequence = CLARNA::FindInteractions::m_sequence;
				//seqWithSeparateChains = CLARNA::FindInteractions::seqWithSeparateChains;
			//}
		//}
	//}
	//else
	//{
		//cout << "12DDSimScore for 3D strucutre only accept directory path which includes all the required files\n";
		//exit(EXIT_FAILURE);
	//}
//}

void CLARNA_all::FindInteractions::readInputFiles(const std::filesystem::path& path)
{
	const std::filesystem::path target_dir { path };
	std::vector<FullClaRNATuple> vfct_pdb, vfct_clarna;
	
	if(std::filesystem::is_directory(path))
	{
		//cout << "finding sequence file in folder and read it\n";
		const std::filesystem::path target_dir { path };
		// iterate all over the target directory to find the sequence file
		for(const auto& dir_entry : std::filesystem::directory_iterator { target_dir })
		{
			const std::filesystem::path tmp_path { dir_entry };
			if(tmp_path.extension() == ".pdb")
			{
				vfct_pdb = CLARNA::FindInteractions::readPDBFile(tmp_path);
				m_sequence = CLARNA::FindInteractions::m_sequence;
				seqWithSeparateChains = CLARNA::FindInteractions::seqWithSeparateChains;
				m_pdbFlag = true;
			}
			if(tmp_path.extension() == m_extension)
			{
				m_extFlag = true;
				m_vStructures_name.emplace_back(tmp_path.string());
				auto vct_clarna = CLARNA::FindInteractions::readInputFile(tmp_path);
				vfct_clarna = CLARNA::FindInteractions::makeReciprocalInteractions(
							  CLARNA::FindInteractions::makeFullVectorOfClaRNATuple(vct_clarna));
							  
				auto vfct =   CLARNA::FindInteractions::addNonInteractedEdges(vfct_pdb, vfct_clarna);
				CLARNA::FindInteractions::vftAll = CLARNA::FindInteractions::finalBasePairCheck(vfct);
				v_vftAll.emplace_back(CLARNA::FindInteractions::vftAll);
				
				CLARNA::FindInteractions::sepCanBasePairs(vct_clarna, vfct_pdb);
				v_vftCans.emplace_back(CLARNA::FindInteractions::vftCans);
				
				CLARNA::FindInteractions::sepNonCanBasePairs(vct_clarna, vfct_pdb);
				v_vftNonCans.emplace_back(CLARNA::FindInteractions::vftNonCans);
				
				CLARNA::FindInteractions::sepWobble(vct_clarna, vfct_pdb);
				v_vftWobbles.emplace_back(CLARNA::FindInteractions::vftWobbles);
				
				CLARNA::FindInteractions::sepStacking(vct_clarna, vfct_pdb);
				v_vftStacks.emplace_back(CLARNA::FindInteractions::vftStacks);
				
				CLARNA::FindInteractions::sepBasePairs(vct_clarna, vfct_pdb);
				v_vftBasePairs.emplace_back(CLARNA::FindInteractions::vftBasePairs);
			}
		}
		if(!m_pdbFlag)
		{
			cout << "PDB file is missed in your dataset\n";
			exit(EXIT_FAILURE);
		}
		
		if(!m_extFlag)
		{
			cout << "There is no file with " << m_extension << " in your dataset\n";
			exit(EXIT_FAILURE);
		}
	}
}


