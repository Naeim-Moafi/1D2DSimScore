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


//=======================================================================================================
BLAST::FindInteractions::FindInteractions()
	: m_numberOfCanBasePairs{0}, m_numberOfNonCanBasePairs{0}, m_numberOfAllBasePairs{0} {}

BLAST::FindInteractions::FindInteractions(bool isWobble_canonical)
	: FindInteractions()
{
	m_isWobble_canonical = isWobble_canonical;
}

//split_string
int BLAST::FindInteractions::split_string(std::string_view string2parse, const char& delim, std::vector<std::string>& output_list) const
{
	int numberOfItems = 0;
	std::string currentToken;
	std::istringstream iss(string2parse.data());
	// get each part of iss which are separated with delim
	while(std::getline(iss, currentToken, delim))
	{
		output_list.push_back(currentToken);
		++numberOfItems;
	}
		
	return numberOfItems;
}

//init_ranges
void BLAST::FindInteractions::init_ranges(std::ifstream& inpFile)
{
	std::string line;
	std::vector<std::string> vStrOut;
	Range range;
	int n = 0;
	inpFile.clear();
	inpFile.seekg(0);
	while(std::getline(inpFile, line))
	{
		if(line == ">ranges")
		{
			while(line.at(0) != '>' || line.size() != 0)
			{
				std::getline(inpFile, line);
				// call getline, we should double-check the conditions once more
				if(line.size() == 0 || line.at(0) == '>')
				{
					break;
				}
				//choosing # character for commenting in the file
				if(line.at(0) == '#')
				{
					continue;
				}
				
				split_string(line, ':', vStrOut);
				range.start = std::stoi(vStrOut[n]);
				if(range.start == 0)
				{
					cerr << "Error: the range should be start from 1\n";
					exit(EXIT_FAILURE);
				}
				range.term = std::stoi(vStrOut[n + 1]);
				m_vRanges.push_back(range);
				n += 2;
			}
		}
	}
}

//init_sequences
void BLAST::FindInteractions::init_sequences(std::ifstream& inpFile)
{
	std::string line;
	std::vector<std::string> vStrOut;
	inpFile.clear();
	inpFile.seekg(0);
	while(std::getline(inpFile, line))
	{
		if(line == ">sequences")
		{
			while(line.at(0) != '>' || line.size() != 0)
			{
				std::getline(inpFile, line);
				// call getline, we should double-check the conditions once more
				if(line.size() == 0 || line.at(0) == '>')
				{
					break;
				}
				//choosing # character for commenting in the file
				if(line.at(0) == '#')
				{
					continue;
				}
					
				m_vSequences.push_back(line);
			}
		}
	}	
}

//init_sequences
void BLAST::FindInteractions::init_structures(std::ifstream& inpFile)
{
	std::string line;
	std::vector<std::string> vStrOut;
	inpFile.clear();
	inpFile.seekg(0);
	while(std::getline(inpFile, line))
	{
		if(line == ">structures")
		{
			while(line.at(0) != '>' || line.size() != 0)
			{
				std::getline(inpFile, line);
				// call getline, we should double-check the conditions once more
				if(line.size() == 0 || line.at(0) == '>')
				{
					break;
				}
				//choosing # character for commenting in the file
				if(line.at(0) == '#')
				{
					continue;
				}
				m_vStructures.push_back(line);
			}
		}
	}
}

// resetIndices
//SSMap BLAST::FindInteractions::resetIndices(Range range, std::pair<int, int> ssPair) const
//{
	//if(ssPair.first > range.start)
	//{
		//if(ssPair.second < range.term)
		//{
			//return 
		//}
	//}
//}

void BLAST::FindInteractions::set_isWobble_canonical(bool isWobble_canonical)
{
	SS::FindInteractions::set_isWobble_canonical(isWobble_canonical);
}

SSMap BLAST::FindInteractions::resetIndices(Range range, SSMap ssMap) const
{
	SSMap ssMap_new;
	for(const auto& ssPair : ssMap)
	{
		if(range.start <= ssPair.first && ssPair.first <= range.term)
		{
			if(ssPair.second <= range.term )
			{
				if( ssPair.second == -1)
				{
					ssMap_new.insert(std::make_pair(ssPair.first - range.start + 1, ssPair.second));
				}
				else
				{
					ssMap_new.insert(std::make_pair(ssPair.first - range.start + 1, ssPair.second - range.start + 1));
				}
			}
			else
			{
				ssMap_new.insert(std::make_pair(ssPair.first - range.start + 1, -3));
			}
		}
		else if(ssPair.first <= range.start)
		{
			if(range.start <= ssPair.second && ssPair.second <= range.term)
			{
				if(ssMap_new.find(ssPair.second) == ssMap_new.end())
				{
					ssMap_new.insert(std::make_pair(ssPair.second - range.start + 1, -2));
				}
			}
		}
	}
	return ssMap_new;
}

std::vector<SSMap> BLAST::FindInteractions::resetIndices(std::vector<Range> vRanges, std::vector<SSMap> vSSMap) const
{
	std::vector<SSMap> vSSMap_new;
	for(size_t i{0}; i < vRanges.size(); ++i)
	{
		vSSMap_new.push_back(resetIndices(vRanges[i], vSSMap[i]));
	}
	return vSSMap_new;
}

//readInputFile
void BLAST::FindInteractions::readInputFile(std::filesystem::path path)
{
	//cout << fmt::format("Start reading file {}", path.filename().string()) << endl;
	//cout << "Start reading file " << path.filename().string() << endl;
	
	//int numberOfLine = 0;
	std::ifstream inpFile(path.string().c_str());
	std::string line;
	
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
	
	init_ranges(inpFile);
	init_sequences(inpFile);
	init_structures(inpFile);
	
	// the size of the m_vRanges, m_vSequences, and m_vStructures must be the same
	if(m_vRanges.size() != m_vSequences.size() || m_vRanges.size() != m_vStructures.size() || m_vSequences.size() != m_vStructures.size())
	{
		cout << "something is wrong" << endl;
		//cout << fmt::format("there are {} ranges, {} sequences, {} structures", m_vRanges.size(),  m_vSequences.size(), m_vStructures.size()) << endl;
		cout << "there are " << m_vRanges.size() << " ranges, " << m_vSequences.size() << " sequences, and " << m_vStructures.size() << " structures" << endl; 
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
	// in this stage all the size are the same
	size_t sizeOfVectors = m_vRanges.size();
	//cout << fmt::format("there are information for {} structures", sizeOfVectors) << endl;
	/*cout << "there are information for " << sizeOfVectors << " structures" << endl;
	for(size_t i{0}; i < sizeOfVectors; ++i)
	{
		cout << m_vRanges[i].start << " " << m_vRanges[i].term  << endl;
		cout << m_vSequences[i] << endl;
		cout << m_vStructures[i] << endl;
	}*/
	inpFile.close();
	
	for( size_t i{0}; i < sizeOfVectors; ++i)
	{
		m_vSSMaps.push_back(SS::FindInteractions::extractInteractions(m_vStructures[i]));
		SS::FindInteractions::mapSS2seq(m_vSequences[i], m_vSSMaps[i]);
		m_vSSMaps_can.push_back(SS::FindInteractions::ssMap_can);
		m_vSSMaps_noncan.push_back(SS::FindInteractions::ssMap_noncan);
		SS::FindInteractions::ssMap_can.clear();
		SS::FindInteractions::ssMap_noncan.clear();
	}
}

std::string BLAST::FindInteractions::structure2Can(const std::string& sequence, const std::string& structure, SSMap ssMap) const
{
	std::string newStruct(structure);
	for(const auto& pair : ssMap)
	{
		if(!SS::FindInteractions::isCanonical(sequence, pair))
		{
			if(pair.second != -1)
			{
				newStruct.replace(pair.first - 1, 1, ".");
				newStruct.replace(pair.second - 1, 1, ".");
			}
		}
	}
	
	return newStruct;
}

void BLAST::FindInteractions::structure2Can() 
{
	for(size_t i{0}; i < m_vRanges.size(); ++i)
	{
		m_vCansStructures.push_back(structure2Can(m_vSequences[i], m_vStructures[i], m_vSSMaps[i]));
	}
}


std::string BLAST::FindInteractions::structure2NonCan(const std::string& sequence, const std::string& structure, SSMap ssMap) const
{
	std::string newStruct(structure);
	for(const auto& pair : ssMap)
	{
		if(SS::FindInteractions::isCanonical(sequence, pair))
		{
			newStruct.replace(pair.first - 1, 1, ".");
			newStruct.replace(pair.second - 1, 1, ".");
		}
	}
	
	return newStruct;
}

void BLAST::FindInteractions::structure2NonCan() 
{
	for(size_t i{0}; i < m_vRanges.size(); ++i)
	{
		m_vNonCansStructures.push_back(structure2NonCan(m_vSequences[i], m_vStructures[i], m_vSSMaps[i]));
	}
}

void BLAST::FindInteractions::showInformation() const
{
	cout << "there are information of " << m_vRanges.size() << "structures" << endl;
	for(size_t i{0}; i < m_vRanges.size(); ++i)
	{
		cout << "Structure #" << i+1 << ":\n";
		cout << "Ranges:\n\t";
	    cout << "Start: " << m_vRanges[i].start << " end: " << m_vRanges[i].term  << endl;
	    cout << "Sequence:\n\t";
	    cout << m_vSequences[i] << endl;
	    cout << "Secondary structure:\n\t";
	    cout << m_vStructures[i] << endl << endl;
	}
}

void BLAST::FindInteractions::showInformationInRange(Range range, SSMap ssMap, int structure_number, const std::string& interaction_type) const
{
	//sequence
	cout << "\t\t";
	for(auto i { range.start }; i <= range.term; ++i)
	{
		cout << m_vSequences[structure_number][i-1];
	}
	cout << endl;
	
	//secondary structure
	cout << "\t\t";
	for(auto i { range.start }; i <= range.term; ++i)
	{
		if(ssMap[(i - range.start) + 1] != -2 && ssMap[(i - range.start) + 1] != -3)
		{
			if(interaction_type == "A" || interaction_type == "a")
			{
				cout << m_vStructures[structure_number][i-1];
			}
			if(interaction_type == "C" || interaction_type == "c")
			{
				cout << m_vCansStructures[structure_number][i-1];
			}
			if(interaction_type == "N" || interaction_type == "n")
			{
				cout << m_vNonCansStructures[structure_number][i-1];
			}
		}
		else
		{
			if(ssMap[(i - range.start) + 1] == -3)
			{
				cout << "x";
			}
			
			if(ssMap[(i - range.start) + 1] == -2)
			{
				cout << "y";
			}
		}
	}
	cout << endl;
}




//#define __MAIN__
#ifdef __MAIN__
using namespace std;
//using namespace SS;
using namespace filesystem;
 
int main()
{	
	BLAST::FindInteractions fiBlast;
	fiBlast.readInputFile("blast_example.txt");
	vector<SSMap> vSSMap = fiBlast.resetIndices(fiBlast.m_vRanges, fiBlast.m_vSSMaps);
	//int n = 0;
	//for(const auto& ssMap : vSSMap)
	//{
		////cout << fmt::format("structure #{}", ++n) << endl;
		//cout << "structure #" << ++n << endl;
 		//for(const auto& ssPair : ssMap)
		//{
			////cout << fmt::format("({:3}, {:3})", ssPair.first, ssPair.second) << endl;
			//cout << ssPair.first << " " << ssPair.second << endl;
		//}
		//cout << endl;
	//}
	
	//fiBlast.structure2Can(fiBlast.m_vSequences[0], fiBlast.m_vStructures[0], fiBlast.m_vSSMaps[0]);
	//getchar();
	fiBlast.structure2Can(fiBlast.m_vSequences[1], fiBlast.m_vStructures[1], fiBlast.m_vSSMaps[1]);
	return 0;
}

#endif //__MAIN__
