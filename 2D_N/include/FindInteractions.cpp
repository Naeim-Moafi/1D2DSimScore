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
									   
									   
									  
//void initClaRNAMatrix(ClaRNAMatrix& cMat, size_t length)
//{
	//cMat.clear();
	//for(size_t i { 0 }; i < length; ++i)
	//{
		//std::vector<ClaRNAMatrix_elem> v_tmp;
		//for(size_t j { 0 }; j < length; ++j)
		//{
			//v_tmp.push_back(std::make_tuple(std::make_tuple<0,0,0,0,0,0>,std::make_tuple<0,0,0,0,0,0>,std::make_tuple<0,0,0,0,0,0>,std::make_tuple<0,0>,std::make_tuple<0,0>));
		//}
		//cMat.emplace_back(v_tmp);
	//}
//}									 

//=======================================================================================================
CLARNA::FindInteractions::FindInteractions()
{
	m_numberOfCanBasePairs = 0;
	m_numberOfNonCanBasePairs = 0;
	m_numberOfAllBasePairs = 0;
	m_numberOfStacks = 0;
	m_numberOfAllInteractions = 0;
	m_number_involved_faces_edges = 5;
	
}

CLARNA::FindInteractions::FindInteractions(bool isWobble_canonical)
	: CLARNA::FindInteractions()
{
	m_isWobble_canonical = isWobble_canonical;	
}


int CLARNA::FindInteractions::get_numberOfCanBasePairs()const
{
	return m_numberOfCanBasePairs/2;
}

int CLARNA::FindInteractions::get_numberOfNonCanBasePairs() const
{
	return m_numberOfNonCanBasePairs/2;
}

int CLARNA::FindInteractions::get_numberOfStacks() const
{
	return m_numberOfStacks/2;
}

int CLARNA::FindInteractions::get_numberOfAllInteractions() const
{
	return m_numberOfAllInteractions;
}

int CLARNA::FindInteractions::get_numberOfAllBasePairs() const
{
	return m_numberOfAllBasePairs/2;
}

bool CLARNA::FindInteractions::get_isWobble_canonical() const
{
	return m_isWobble_canonical;
}

void CLARNA::FindInteractions::set_isWobble_canonical(bool isWobble_canonical)
{
	m_isWobble_canonical = isWobble_canonical;
}

int CLARNA::FindInteractions::get_number_involved_faces_edges() const
{
	return m_number_involved_faces_edges;
}

void CLARNA::FindInteractions::set_number_involved_faces_edges(int number_involved_faces_edges)
{
	m_number_involved_faces_edges = number_involved_faces_edges;
}

ClaRNATuple CLARNA::FindInteractions::extractInteractions(std::string clarnaLine)
{
	ClaRNATuple ct;
	std::istringstream iss(clarnaLine);
	std::string word;
	std::vector<std::string> vStrData; // all the information of each line of the clarna out would be saved in vector 
									   // and then will be parsed into ClaRNATuple.
	while(iss >> word)
	{
		vStrData.push_back(word);
	}
	
	
	ct = std::make_tuple(static_cast<std::string>(vStrData[0]), stoi(vStrData[1]), static_cast<std::string>(vStrData[2]), static_cast<std::string>(vStrData[3]), stoi(vStrData[4]), static_cast<std::string>(vStrData[5]), static_cast<std::string>(vStrData[6]), stod(vStrData[7]));
	return ct;
}

std::vector<ClaRNATuple> CLARNA::FindInteractions::readInputFile(std::filesystem::path path)
{
	std::vector<ClaRNATuple> vct;
	std::string currLine;
	
	//cout << fmt::format("start reading file {}", path.filename().string()) << endl;
	//cout << "start reading file " <<  path.filename().string() << endl;
	
	
	// open the file in path
	std::ifstream inpFile(path.string().c_str());
	if(inpFile.fail())
	{
		if(!path.parent_path().empty())
		{
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), path.relative_path().string()) << endl;
			cout << "Please check the file '" << path.filename().string() << "' in the followng directory\n'" << path.relative_path().string() << "'" << endl;
		}
		else
		{
			std::filesystem::path p = std::filesystem::current_path();
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), p.relative_path().string()) << endl;
			cout << "Please check the file '" << p.filename().string() << "' in the followng directory\n'" << p.relative_path().string() << "'" << endl;
		}
		exit(EXIT_FAILURE);
	}
	
	// check the input file line by line and extract the interactions
	while(!inpFile.fail())
	{
		std::getline(inpFile, currLine);
		//in clarna output the last line is empty and without any information
		if(currLine.size() == 0)
		{
			//it is the last line;
			break;
		}
		vct.push_back(extractInteractions(currLine));
		++m_numberOfAllInteractions;
		
	}
	
	
	inpFile.close();
	
	return vct;
}

// make full ClaRNATuple
FullClaRNATuple CLARNA::FindInteractions::makeFullClaRNATuple(const ClaRNATuple& ct)
{
	std::string type; // c(cis) or t(trans) for basepirs and empty for stacking
	FullClaRNATuple fct;
	auto[firstChainID, firstNuclNumber, firstNuclName, secondChainID, secondNuclNumber, secondNuclName, typeOfInteraction, weight] = ct;
	std::string firstEdge = typeOfInteraction.substr(0,1);
	std::string secondEdge = typeOfInteraction.substr(1,1);
	
	//if(isBaseBackbone(ct))
	//{
		//stdExt::tuplePrint(ct);
		//if(isBasePhosphate(ct))
		//{
			//cout << "BP" <<endl;
			
		//}
		//if(isBaseRibose(ct))
		//{
			//cout << "BR" << endl;
		//}
		//getchar();
	//}
	
	
	if(isStacking(ct))
	{
		type = ""; // there is no cis or trans for stackings, >> >< <> <<
	}
	else
	{
		type = typeOfInteraction.substr(3,1); // c or t;
	}
	
	fct = std::make_tuple(firstChainID, firstNuclNumber, firstNuclName, firstEdge, secondChainID, secondNuclNumber, secondNuclName, secondEdge, type, weight);
	return fct;
}

// makeFullVectorOfClaRNATuple
std::vector<FullClaRNATuple> CLARNA::FindInteractions::makeFullVectorOfClaRNATuple(const std::vector<ClaRNATuple>& vct)
{
	std::vector<FullClaRNATuple> vfct;
	// iterate through the vct and make a FullClaRNATuple and push it back to the vfct
	for(const auto& ct : vct)
	{
		vfct.push_back(makeFullClaRNATuple(ct));
	}
	
	return vfct;
}


//stackingReciprocal
std::string CLARNA::FindInteractions::stackingReciprocal(const std::string nuclEdge) const
{
	if(nuclEdge == ">")
	{
		return "<";
	}
	
	return ">";
}

//makeReciprocalInteractions
std::vector<FullClaRNATuple> CLARNA::FindInteractions::makeReciprocalInteractions(const  std::vector<FullClaRNATuple>& vfct) const
{
	std::vector<FullClaRNATuple> vfct_new;
	for(const auto& fct : vfct)
	{
		auto [firstChainID, firstNuclNumber, firstNuclName, firstEdge, secondChainID, secondNuclNumber, secondNuclName, secondEdge, type, weight] = fct;
		vfct_new.push_back(fct);
		if(isStacking(firstEdge))
		{
			vfct_new.push_back(std::make_tuple(secondChainID, secondNuclNumber, secondNuclName, stackingReciprocal(secondEdge), firstChainID, firstNuclNumber, firstNuclName, stackingReciprocal(firstEdge), type, weight));
		}
		else
		{
			vfct_new.push_back(std::make_tuple(secondChainID, secondNuclNumber, secondNuclName, secondEdge, firstChainID, firstNuclNumber, firstNuclName, firstEdge, type, weight));
		}
	}
	
	
	return vfct_new;
}

//readPDBFile
std::vector<FullClaRNATuple> CLARNA::FindInteractions::readPDBFile(std::filesystem::path path)
{
	m_vNucleotides.clear();
	ReadPDB inputPDB;
	std::vector<std::string> vSeqOfChains;
	std::ostringstream ossSeq;
	std::ostringstream oSeq;
	//cout << fmt::format("start reading {} file", path.filename().string()) << endl;
	////cout << "start reading  " <<  path.filename().string() << endl;
	inputPDB.read_file(path);
	int numberOfNucleotides = inputPDB.get_n_residues();
	int numberOfNucleotidesInCurrChain;
	int numberOfChains = inputPDB.get_n_chains();
	std::vector<Range> chains(numberOfChains);
	Atom nucleotides[numberOfNucleotides];
	// we need to reserved capacity in the amount of the number of nucleotide for vNucleotides
	m_vNucleotides.reserve(numberOfNucleotides);
	int numberOfErrors{0};
	// check all nucleotide in each chain and extract information 
	for(int i{0}; i < numberOfChains; ++i)
	{
		numberOfNucleotidesInCurrChain = inputPDB.get_n_residues_in_chain(i);
		inputPDB.select_chain_for_reading(i);
		for(int j{0}; j < numberOfNucleotidesInCurrChain; ++j)
		{
			// using atom C4' for extracting informations
			int isReadC4_OK = inputPDB.get_atom_from_selected_chain(nucleotides[chains[i].start + j], j, "C4'|C4*");
			// check if the pdb file has gaps or not
			if(!isReadC4_OK)
			{ // if yes report the gaps
				//cout << fmt::format("the atom C4' in {}th nucleotide in chain {} is missed.", j, possibleChainID[i]) << endl;
				cout << "the atom C4' in " << j << "th nucleotide in chain " << possibleChainID[i] << " is missed."  << endl;
				cout << "Skipping to initializing this nucleotid." << endl;
				++numberOfErrors;
				continue;
			}
			ossSeq << nucleotides[chains[i].start + j].nuclName;
			m_sequence += nucleotides[chains[i].start + j].nuclName;
			m_vNucleotides.push_back(nucleotides[chains[i].start + j]);			
		}
		
		vSeqOfChains.push_back(ossSeq.str());
		ossSeq.str("");
		ossSeq.clear();
	};
	
	if(numberOfErrors > 0)
	{
		if(!path.parent_path().empty())
		{
			cout << "the program found " << numberOfErrors << " gap(s) in " << path.filename().string() << " in the following directory\n" << path.relative_path().string() << endl;
		}
		else
		{
			std::filesystem::path p = std::filesystem::current_path();
			//cout << fmt::format("the program found {} gap(s) in {} in the following directory\n{}", numberOfErrors, path.filename().string(), p.relative_path().string()) << endl;
			cout << "the program found " << numberOfErrors << " gap(s) in " << p.filename().string() << " in the following directory\n" << p.relative_path().string() << endl;
			
		}
		exit(EXIT_FAILURE);
	}
	std::copy(std::cbegin(vSeqOfChains), std::cend(vSeqOfChains), std::ostream_iterator<std::string>(oSeq, " "));
	seqWithSeparateChains = oSeq.str();
	std::vector<FullClaRNATuple> vfct;
	switch(m_number_involved_faces_edges)
	{
		case 1 :  vfct = makeFullyNonInteractedVector1(); break;
		case 2 :  vfct = makeFullyNonInteractedVector2(); break;
		case 3 :  vfct = makeFullyNonInteractedVector3(); break;
		case 5 :  vfct = makeFullyNonInteractedVector5(); break;
		default : cout << "Incorrect number of faces and edges are requested\n"<< endl;
	}
	
	
	return vfct;
}

// makeFullyNonInteractedVector
std::vector<FullClaRNATuple> CLARNA::FindInteractions::makeFullyNonInteractedVector5()
{
	
	std::vector<FullClaRNATuple> vfct;
	
	for(const auto& nucleotide : m_vNucleotides)
	{
		std::string sChainID(nucleotide.chainId);
		std::string sNuclName(nucleotide.nuclName);
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "W", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "S", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "H", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, ">", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "<", "_", -1, "_", "_", "_", 0.0));
	}
	
	return vfct;
}

std::vector<FullClaRNATuple> CLARNA::FindInteractions::makeFullyNonInteractedVector1()
{
	
	std::vector<FullClaRNATuple> vfct;
	
	for(const auto& nucleotide : m_vNucleotides)
	{
		std::string sChainID(nucleotide.chainId);
		std::string sNuclName(nucleotide.nuclName);
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "W", "_", -1, "_", "_", "_", 0.0));
	}
	
	return vfct;
}

std::vector<FullClaRNATuple> CLARNA::FindInteractions::makeFullyNonInteractedVector3()
{
	
	std::vector<FullClaRNATuple> vfct;
	
	for(const auto& nucleotide : m_vNucleotides)
	{
		std::string sChainID(nucleotide.chainId);
		std::string sNuclName(nucleotide.nuclName);
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "W", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "H", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "S", "_", -1, "_", "_", "_", 0.0));
	}
	
	return vfct;
}


std::vector<FullClaRNATuple> CLARNA::FindInteractions::makeFullyNonInteractedVector2()
{
	
	std::vector<FullClaRNATuple> vfct;
	
	for(const auto& nucleotide : m_vNucleotides)
	{
		std::string sChainID(nucleotide.chainId);
		std::string sNuclName(nucleotide.nuclName);
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, ">", "_", -1, "_", "_", "_", 0.0));
		vfct.push_back(std::make_tuple(sChainID, nucleotide.nuclNumber, sNuclName, "<", "_", -1, "_", "_", "_", 0.0));
	}
	
	return vfct;
}


bool CLARNA::FindInteractions::isCanonical(const ClaRNATuple& ct) const 
{
	//extract information of the ct
	std::string firstNuclName =  std::get<2>(ct);
	std::string secondNuclName = std::get<5>(ct);
	std::string typeOfInteraction = std::get<6>(ct);	
	
	// all of the canonical interaction are WW_cis
	if(typeOfInteraction == "WW_cis")
	{
		// no we can check the base pair
		if( (firstNuclName == "A" && secondNuclName == "U") || (firstNuclName == "U" && secondNuclName == "A") )
		{
			return true;
		}
		
		if( (firstNuclName == "C" && secondNuclName == "G") || (firstNuclName == "G" && secondNuclName == "C") )
		{
			return true;
		}
		
		// it deponds on user if wobble is canonical or not
		if(get_isWobble_canonical())
		{
			if( (firstNuclName == "G" && secondNuclName == "U") || (firstNuclName == "U" && secondNuclName == "G"))
			{
				return true;
			}
		}
	}
	
	return false;
}

bool CLARNA::FindInteractions::isCanonical(const FullClaRNATuple& fct) const
{
	
	//extract required information from the ct
	std::string firstNuclName = std::get<2>(fct);
	std::string firstEdge = std::get<3>(fct);
	std::string secondNuclName = std::get<6>(fct);
	std::string secondEdge = std::get<7>(fct);
	std::string type = std::get<8>(fct);
	
	// canonical interactions are AU, UA, CG, GC with WW_cis type
	// checking the edges which are involved in the interactions
	if(firstEdge == "W" && secondEdge == "W" && type == "c")
	{
		// if edges are correct for canonical interactions
		// we need to check the nucleotides which are invovleved
		if((firstNuclName == "A" && secondNuclName == "U") || (firstNuclName == "U" && secondNuclName == "A")) 
		{
			return true;
		}
		
		if((firstNuclName == "C" && secondNuclName == "G") || (firstNuclName == "G" && secondNuclName == "C")) 
		{
			return true;
		}
		
		// checking user request
		if(get_isWobble_canonical())
		{
			if((firstNuclName == "U" && secondNuclName == "G") || (firstNuclName == "G" && secondNuclName == "U"))
			{
				return true;
			}
		}
	}
	
	return false;
}


bool CLARNA::FindInteractions::isWobble(const ClaRNATuple& ct) const 
{
	//extract information of the ct
	std::string firstNuclName =  std::get<2>(ct);
	std::string secondNuclName = std::get<5>(ct);
	std::string typeOfInteraction = std::get<6>(ct);	
	
	// all of the canonical interaction are WW_cis
	if(typeOfInteraction == "WW_cis")
	{
		if( (firstNuclName == "G" && secondNuclName == "U") || (firstNuclName == "U" && secondNuclName == "G"))
		{
			return true;
		}
	
	}
	
	return false;
}


bool CLARNA::FindInteractions::isWobble(const FullClaRNATuple& fct) const
{
	
	//extract required information from the ct
	std::string firstNuclName = std::get<2>(fct);
	std::string firstEdge = std::get<3>(fct);
	std::string secondNuclName = std::get<6>(fct);
	std::string secondEdge = std::get<7>(fct);
	std::string type = std::get<8>(fct);
	
	// canonical interactions are AU, UA, CG, GC with WW_cis type
	// checking the edges which are involved in the interactions
	if(firstEdge == "W" && secondEdge == "W" && type == "c")
	{
		if((firstNuclName == "U" && secondNuclName == "G") || (firstNuclName == "G" && secondNuclName == "U"))
		{
			return true;
		}

	}
	
	return false;
}

// isStacking
bool CLARNA::FindInteractions::isStacking(const ClaRNATuple& ct) const
{
	//extract required information from the ct
	std::string typeOfInteraction = std::get<6>(ct);
	if(	   typeOfInteraction == ">>"
		|| typeOfInteraction == "<<"
		|| typeOfInteraction == "><"
		|| typeOfInteraction == "<>")
	{
		return true;
	}
	
	return false;
}

bool CLARNA::FindInteractions::isStacking(const std::string nuclEdge) const
{
	if(nuclEdge == ">" or nuclEdge == "<")
	{
		return true;
	}
	
	return false;
}

//isBaseBackbone
bool CLARNA::FindInteractions::isBasePhosphate(const ClaRNATuple& ct) const
{
	// extract required information from the ct
	std::string  typeOfInteraction = std::get<6>(ct);
	auto itPh = std::find_if(std::begin(vBPh), std::end(vBPh), [typeOfInteraction](const std::string& str){return str == typeOfInteraction;});
	
	if(itPh != std::end(vBPh))
	{
		return true;
	}
	
	return false;
}

bool CLARNA::FindInteractions::isBaseRibose(const ClaRNATuple& ct) const
{
	// extract required information from the ct
	std::string  typeOfInteraction = std::get<6>(ct);
	auto itR = std::find_if(std::begin(vBR), std::end(vBR), [typeOfInteraction](const std::string& str){return str == typeOfInteraction;});
	
	if(itR != std::end(vBR))
	{
		return true;
	}
	
	return false;
	
}

bool CLARNA::FindInteractions::isBaseBackbone(const ClaRNATuple& ct) const
{
	if(isBasePhosphate(ct) || isBaseRibose(ct))
	{
		return true;
	}
	
	return false;
}

// addNonInteractedEdges
std::vector<FullClaRNATuple> CLARNA::FindInteractions::addNonInteractedEdges(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_fromClarna) const
{
	if(vfct_fromClarna.size() == 0)
	{
		return vfct_fromPDB;
	}

	std::vector<FullClaRNATuple> vfct_nonInteractedEdges;
	// insided clarna outputs, it is possible one edge involved in several interactions
	// but 1DSimScore is interested in the one with highest weight
	// the following lambda function does that
	auto findTheBestAccuracy { [](auto vIts){
											FullClaRNATuple fct;
											double maxWeight = 0;
											for(const auto& it : vIts)
											{
												double weight = std::get<9>(*it);
												if(weight > maxWeight)
												{
													maxWeight = weight;
													fct = *it;
												}
											}
											return fct;
										  }
						   };
				   

// the vfct_fromPDB contains all of the nucleotides inside pdb.
	// then we iterate through the vfct_fromPDB 
	// to see if clarna has detected any type of ineractions for the correspoding nucleotide or not
	for(const auto& fct_fromPDB : vfct_fromPDB)
	{
		// check chainID, nuclNumber, nuclName, Edge
		auto nuclEdgeMatches = stdExt::find_all_if(std::begin(vfct_fromClarna), std::end(vfct_fromClarna), [&fct_fromPDB](auto fct_fromClarna){
																										      return (   std::get<0>(fct_fromClarna) == std::get<0>(fct_fromPDB)
																													  && std::get<1>(fct_fromClarna) == std::get<1>(fct_fromPDB)
																													  && std::get<2>(fct_fromClarna) == std::get<2>(fct_fromPDB)
																													  && std::get<3>(fct_fromClarna) == std::get<3>(fct_fromPDB)
																													 );
																											 }
							                      );
		
		if(nuclEdgeMatches.size() > 0)
		{
			vfct_nonInteractedEdges.push_back(findTheBestAccuracy(nuclEdgeMatches));
		}
		else
		{
			vfct_nonInteractedEdges.push_back(fct_fromPDB);
		}
	}
	
	return vfct_nonInteractedEdges;
}

void CLARNA::FindInteractions::addNonInteractedEdges(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath)
{
	if(m_number_involved_faces_edges < FIVE_EDGES)
	{
		cout << "ERROR: Incorrect number of involved faces and edges for all interactions\n";
		cout << "you can use 5\n";
		exit(EXIT_FAILURE);
	}
	
	std::vector<FullClaRNATuple> vfct_clarna = makeReciprocalInteractions(makeFullVectorOfClaRNATuple(readInputFile(clarnaPath)));
	std::vector<FullClaRNATuple> vfct_pdb = readPDBFile(pdbPath);
	std::sort(std::begin(vfct_clarna), std::end(vfct_clarna));
	std::sort(std::begin(vfct_pdb), std::end(vfct_pdb));
	std::vector<FullClaRNATuple> vfct = addNonInteractedEdges(vfct_pdb, vfct_clarna);
	vftAll = finalBasePairCheck(vfct);
	checkClaRNAwithPDB(vftAll, vfct_pdb);
}

std::vector<FullClaRNATuple> CLARNA::FindInteractions::addNonInteractedEdgesForCans(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_fromClarna) const
{
	if(vfct_fromClarna.size() == 0)
	{
		return vfct_fromPDB;
	}
	std::vector<FullClaRNATuple> vfct_nonInteractedEdges;
	// insided clarna outputs, it is possible one edge involved in several interactions
	// but 1DSimScore is interested in the one with highest weight
	// the following lambda function does that
	auto findTheBestAccuracy { [](auto vIts){
											FullClaRNATuple fct;
											double maxWeight = 0;
											for(const auto& it : vIts)
											{
												double weight = std::get<9>(*it);
												if(weight > maxWeight)
												{
													maxWeight = weight;
													fct = *it;
												}
											}
											return fct;
										  }
						   };
				   

// the vfct_fromPDB contains all of the nucleotides inside pdb.
	// then we iterate through the vfct_fromPDB 
	// to see if clarna has detected any type of ineractions for the correspoding nucleotide or not
	for(const auto& fct_fromPDB : vfct_fromPDB)
	{
		// check chainID, nuclNumber, nuclName, Edge
		auto nuclEdgeMatches = stdExt::find_all_if(std::begin(vfct_fromClarna), std::end(vfct_fromClarna), [&fct_fromPDB](auto fct_fromClarna){
																										      return (   std::get<0>(fct_fromClarna) == std::get<0>(fct_fromPDB)
																													  && std::get<1>(fct_fromClarna) == std::get<1>(fct_fromPDB)
																													  && std::get<2>(fct_fromClarna) == std::get<2>(fct_fromPDB)
																													  && std::get<3>(fct_fromClarna) == std::get<3>(fct_fromPDB)
																													 );
																											 }
							                      );
		
		if(nuclEdgeMatches.size() > 0)
		{
			auto theBestAccuracy = findTheBestAccuracy(nuclEdgeMatches);
			if(isCanonical(theBestAccuracy))
			{
				vfct_nonInteractedEdges.emplace_back(findTheBestAccuracy(nuclEdgeMatches));
			}
			else
			{
				vfct_nonInteractedEdges.emplace_back(fct_fromPDB);
			}
		}
		else
		{
			vfct_nonInteractedEdges.push_back(fct_fromPDB);
		}
	}
	

	
	return vfct_nonInteractedEdges;
}

void CLARNA::FindInteractions::addNonInteractedEdgesForCans(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath)
{		
	if(m_number_involved_faces_edges == TWO_EDGES)
	{
		cout << "ERROR: Incorrect number of involved faces and edges for canonicals\n";
		cout << "you can use 1 or 3 or 5\n";
		exit(EXIT_FAILURE);
	}
	
	std::vector<FullClaRNATuple> vfct_clarna = makeReciprocalInteractions(makeFullVectorOfClaRNATuple(readInputFile(clarnaPath)));
	std::vector<FullClaRNATuple> vfct_pdb = readPDBFile(pdbPath);
	std::sort(std::begin(vfct_clarna), std::end(vfct_clarna));
	std::sort(std::begin(vfct_pdb), std::end(vfct_pdb));
	std::vector<FullClaRNATuple> vfct = addNonInteractedEdgesForCans(vfct_pdb, vfct_clarna);
	vftCans = finalBasePairCheck(vfct);
	checkClaRNAwithPDB(vftCans, vfct_pdb);
}

std::vector<FullClaRNATuple> CLARNA::FindInteractions::addNonInteractedEdgesForWobble(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_fromClarna) const
{
	if(vfct_fromClarna.size() == 0)
	{
		return vfct_fromPDB;
	}
	std::vector<FullClaRNATuple> vfct_nonInteractedEdges;
	// insided clarna outputs, it is possible one edge involved in several interactions
	// but 1DSimScore is interested in the one with highest weight
	// the following lambda function does that
	auto findTheBestAccuracy { [](auto vIts){
											FullClaRNATuple fct;
											double maxWeight = 0;
											for(const auto& it : vIts)
											{
												double weight = std::get<9>(*it);
												if(weight > maxWeight)
												{
													maxWeight = weight;
													fct = *it;
												}
											}
											return fct;
										  }
						   };
				   

// the vfct_fromPDB contains all of the nucleotides inside pdb.
	// then we iterate through the vfct_fromPDB 
	// to see if clarna has detected any type of ineractions for the correspoding nucleotide or not
	for(const auto& fct_fromPDB : vfct_fromPDB)
	{
		// check chainID, nuclNumber, nuclName, Edge
		auto nuclEdgeMatches = stdExt::find_all_if(std::begin(vfct_fromClarna), std::end(vfct_fromClarna), [&fct_fromPDB](auto fct_fromClarna){
																										      return (   std::get<0>(fct_fromClarna) == std::get<0>(fct_fromPDB)
																													  && std::get<1>(fct_fromClarna) == std::get<1>(fct_fromPDB)
																													  && std::get<2>(fct_fromClarna) == std::get<2>(fct_fromPDB)
																													  && std::get<3>(fct_fromClarna) == std::get<3>(fct_fromPDB)
																													 );
																											 }
							                      );
		
		if(nuclEdgeMatches.size() > 0)
		{
			auto theBestAccuracy = findTheBestAccuracy(nuclEdgeMatches);
			if(isWobble(theBestAccuracy))
			{
				vfct_nonInteractedEdges.emplace_back(findTheBestAccuracy(nuclEdgeMatches));
			}
			else
			{
				vfct_nonInteractedEdges.emplace_back(fct_fromPDB);
			}
		}
		else
		{
			vfct_nonInteractedEdges.push_back(fct_fromPDB);
		}
	}
	

	
	return vfct_nonInteractedEdges;
}

void CLARNA::FindInteractions::addNonInteractedEdgesForWobble(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath)
{		
	if(m_number_involved_faces_edges == TWO_EDGES)
	{
		cout << "ERROR: Incorrect number of involved faces and edges for wobbles\n";
		cout << "you can use 1 or 3 or 5\n";
		exit(EXIT_FAILURE);
	}
	
	std::vector<FullClaRNATuple> vfct_clarna = makeReciprocalInteractions(makeFullVectorOfClaRNATuple(readInputFile(clarnaPath)));
	std::vector<FullClaRNATuple> vfct_pdb = readPDBFile(pdbPath);
	std::sort(std::begin(vfct_clarna), std::end(vfct_clarna));
	std::sort(std::begin(vfct_pdb), std::end(vfct_pdb));
	std::vector<FullClaRNATuple> vfct = addNonInteractedEdgesForWobble(vfct_pdb, vfct_clarna);
	vftWobbles = finalBasePairCheck(vfct);
	checkClaRNAwithPDB(vftWobbles, vfct_pdb);
}

std::vector<FullClaRNATuple> CLARNA::FindInteractions::addNonInteractedEdgesForNonCans(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_fromClarna) const
{
	if(vfct_fromClarna.size() == 0)
	{
		return vfct_fromPDB;
	}
	std::vector<FullClaRNATuple> vfct_nonInteractedEdges;
	// insided clarna outputs, it is possible one edge involved in several interactions
	// but 1DSimScore is interested in the one with highest weight
	// the following lambda function does that
	auto findTheBestAccuracy { [](auto vIts){
											FullClaRNATuple fct;
											double maxWeight = 0;
											for(const auto& it : vIts)
											{
												double weight = std::get<9>(*it);
												if(weight > maxWeight)
												{
													maxWeight = weight;
													fct = *it;
												}
											}
											return fct;
										  }
						   };
				   

// the vfct_fromPDB contains all of the nucleotides inside pdb.
	// then we iterate through the vfct_fromPDB 
	// to see if clarna has detected any type of ineractions for the correspoding nucleotide or not
	for(const auto& fct_fromPDB : vfct_fromPDB)
	{
		// check chainID, nuclNumber, nuclName, Edge
		auto nuclEdgeMatches = stdExt::find_all_if(std::begin(vfct_fromClarna), std::end(vfct_fromClarna), [&fct_fromPDB](auto fct_fromClarna){
																										      return (   std::get<0>(fct_fromClarna) == std::get<0>(fct_fromPDB)
																													  && std::get<1>(fct_fromClarna) == std::get<1>(fct_fromPDB)
																													  && std::get<2>(fct_fromClarna) == std::get<2>(fct_fromPDB)
																													  && std::get<3>(fct_fromClarna) == std::get<3>(fct_fromPDB)
																													 );
																											 }
							                      );
		
		if(nuclEdgeMatches.size() > 0)
		{
			auto theBestAccuracy = findTheBestAccuracy(nuclEdgeMatches);
			if(!isCanonical(theBestAccuracy) && !isStacking(std::get<3>(theBestAccuracy)))
			{
				vfct_nonInteractedEdges.emplace_back(theBestAccuracy);
			}
			else
			{
				vfct_nonInteractedEdges.emplace_back(fct_fromPDB);
			}
		}
		else
		{
			vfct_nonInteractedEdges.push_back(fct_fromPDB);
		}
	}
	
	return vfct_nonInteractedEdges;
}

void CLARNA::FindInteractions::addNonInteractedEdgesForNonCans(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath)
{	
	if(m_number_involved_faces_edges == TWO_EDGES || m_number_involved_faces_edges == ONE_EDGE)
	{
		cout << "ERROR: Incorrect number of involved faces and edges for noncanonical base pairs\n";
		cout << "you can use 3 or 5\n";
		exit(EXIT_FAILURE);
	}
	
	std::vector<FullClaRNATuple> vfct_clarna = makeReciprocalInteractions(makeFullVectorOfClaRNATuple(readInputFile(clarnaPath)));
	std::vector<FullClaRNATuple> vfct_pdb = readPDBFile(pdbPath);
	std::sort(std::begin(vfct_clarna), std::end(vfct_clarna));
	std::sort(std::begin(vfct_pdb), std::end(vfct_pdb));
	std::vector<FullClaRNATuple> vfct = addNonInteractedEdgesForNonCans(vfct_pdb, vfct_clarna);
	vftNonCans = finalBasePairCheck(vfct);
	checkClaRNAwithPDB(vftNonCans, vfct_pdb);
}


std::vector<FullClaRNATuple> CLARNA::FindInteractions::addNonInteractedEdgesForStacks(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_fromClarna) const
{
	
	if(vfct_fromClarna.size() == 0)
	{
		return vfct_fromPDB;
	}
	std::vector<FullClaRNATuple> vfct_nonInteractedEdges;
	// insided clarna outputs, it is possible one edge involved in several interactions
	// but 1DSimScore is interested in the one with highest weight
	// the following lambda function does that
	auto findTheBestAccuracy { [](auto vIts){
											FullClaRNATuple fct;
											double maxWeight = 0;
											for(const auto& it : vIts)
											{
												double weight = std::get<9>(*it);
												if(weight > maxWeight)
												{
													maxWeight = weight;
													fct = *it;
												}
											}
											return fct;
										  }
						   };
				   

// the vfct_fromPDB contains all of the nucleotides inside pdb.
	// then we iterate through the vfct_fromPDB 
	// to see if clarna has detected any type of ineractions for the correspoding nucleotide or not
	for(const auto& fct_fromPDB : vfct_fromPDB)
	{
		// check chainID, nuclNumber, nuclName, Edge
		auto nuclEdgeMatches = stdExt::find_all_if(std::begin(vfct_fromClarna), std::end(vfct_fromClarna), [&fct_fromPDB](auto fct_fromClarna){
																										      return (   std::get<0>(fct_fromClarna) == std::get<0>(fct_fromPDB)
																													  && std::get<1>(fct_fromClarna) == std::get<1>(fct_fromPDB)
																													  && std::get<2>(fct_fromClarna) == std::get<2>(fct_fromPDB)
																													  && std::get<3>(fct_fromClarna) == std::get<3>(fct_fromPDB)
																													 );
																											 }
							                      );
		
		if(nuclEdgeMatches.size() > 0)
		{
			auto theBestAccuracy = findTheBestAccuracy(nuclEdgeMatches);
			if(isStacking(std::get<3>(theBestAccuracy)))
			{
				vfct_nonInteractedEdges.emplace_back(theBestAccuracy);
			}
			else
			{
				vfct_nonInteractedEdges.emplace_back(fct_fromPDB);
			}
		}
		else
		{
			vfct_nonInteractedEdges.push_back(fct_fromPDB);
		}
	}
	
	return vfct_nonInteractedEdges;
}

void CLARNA::FindInteractions::addNonInteractedEdgesForStacks(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath)
{	
	if(m_number_involved_faces_edges != TWO_EDGES || m_number_involved_faces_edges != FIVE_EDGES)
	{
		cout << "ERROR: Incorrect number of involved faces and edges for stacking\n";
		cout << "you can use 2 or 5\n";
		exit(EXIT_FAILURE);
	}
	std::vector<FullClaRNATuple> vfct_clarna = makeReciprocalInteractions(makeFullVectorOfClaRNATuple(readInputFile(clarnaPath)));
	std::vector<FullClaRNATuple> vfct_pdb = readPDBFile(pdbPath);
	std::sort(std::begin(vfct_clarna), std::end(vfct_clarna));
	std::sort(std::begin(vfct_pdb), std::end(vfct_pdb));
	std::vector<FullClaRNATuple> vfct = addNonInteractedEdgesForStacks(vfct_pdb, vfct_clarna);
	vftStacks = finalBasePairCheck(vfct);
	checkClaRNAwithPDB(vftStacks, vfct_pdb);
}


std::vector<FullClaRNATuple> CLARNA::FindInteractions::addNonInteractedEdgesForBasePairs(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_fromClarna) const
{
	if(vfct_fromClarna.size() == 0)
	{
		return vfct_fromPDB;
	}
	std::vector<FullClaRNATuple> vfct_nonInteractedEdges;
	// insided clarna outputs, it is possible one edge involved in several interactions
	// but 1DSimScore is interested in the one with highest weight
	// the following lambda function does that
	auto findTheBestAccuracy { [](auto vIts){
											FullClaRNATuple fct;
											double maxWeight = 0;
											for(const auto& it : vIts)
											{
												double weight = std::get<9>(*it);
												if(weight > maxWeight)
												{
													maxWeight = weight;
													fct = *it;
												}
											}
											return fct;
										  }
						   };
				   

// the vfct_fromPDB contains all of the nucleotides inside pdb.
	// then we iterate through the vfct_fromPDB 
	// to see if clarna has detected any type of ineractions for the correspoding nucleotide or not
	for(const auto& fct_fromPDB : vfct_fromPDB)
	{
		// check chainID, nuclNumber, nuclName, Edge
		auto nuclEdgeMatches = stdExt::find_all_if(std::begin(vfct_fromClarna), std::end(vfct_fromClarna), [&fct_fromPDB](auto fct_fromClarna){
																										      return (   std::get<0>(fct_fromClarna) == std::get<0>(fct_fromPDB)
																													  && std::get<1>(fct_fromClarna) == std::get<1>(fct_fromPDB)
																													  && std::get<2>(fct_fromClarna) == std::get<2>(fct_fromPDB)
																													  && std::get<3>(fct_fromClarna) == std::get<3>(fct_fromPDB)
																													 );
																											 }
							                      );
		
		if(nuclEdgeMatches.size() > 0)
		{
			auto theBestAccuracy = findTheBestAccuracy(nuclEdgeMatches);
			if(!isStacking(std::get<3>(theBestAccuracy)))
			{
				vfct_nonInteractedEdges.emplace_back(theBestAccuracy);
			}
			else
			{
				vfct_nonInteractedEdges.emplace_back(fct_fromPDB);
			}
		}
		else
		{
			vfct_nonInteractedEdges.push_back(fct_fromPDB);
		}
	}
	
	return vfct_nonInteractedEdges;
}

void CLARNA::FindInteractions::addNonInteractedEdgesForBasePairs(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath)
{
	if(m_number_involved_faces_edges == TWO_EDGES || m_number_involved_faces_edges == ONE_EDGE)
	{
		cout << "ERROR: Incorrect number of involved faces and edges for all base pairs\n";
		cout << "you can use 3 or 5\n";
		exit(EXIT_FAILURE);
	}
	
	std::vector<FullClaRNATuple> vfct_clarna = makeReciprocalInteractions(makeFullVectorOfClaRNATuple(readInputFile(clarnaPath)));
	std::vector<FullClaRNATuple> vfct_pdb = readPDBFile(pdbPath);
	std::sort(std::begin(vfct_clarna), std::end(vfct_clarna));
	std::sort(std::begin(vfct_pdb), std::end(vfct_pdb));
	std::vector<FullClaRNATuple> vfct = addNonInteractedEdgesForBasePairs(vfct_pdb, vfct_clarna);
	vftBasePairs = finalBasePairCheck(vfct);
	checkClaRNAwithPDB(vftBasePairs, vfct_pdb);
}


void CLARNA::FindInteractions::sepCanBasePairs(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB)
{
	
	std::vector<FullClaRNATuple> vfct;
	// iterate all over the vtInteractions
	for(const auto& ct : vct)
	{
		if(isCanonical(ct))
		{
			vtCans.push_back(ct);
			vfct.push_back(makeFullClaRNATuple(ct));
			//++numberOfCanBasePairs;
		}
	}
	
	auto vfct2 = makeReciprocalInteractions(vfct);
	vftCans = finalBasePairCheck(vfct2);
	vftCans = addNonInteractedEdgesForCans(vfct_formPDB, vfct2);
	checkClaRNAwithPDB(vftCans, vfct_formPDB);
}

void CLARNA::FindInteractions::sepWobble(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB)
{
	
	std::vector<FullClaRNATuple> vfct;
	// iterate all over the vtInteractions
	for(const auto& ct : vct)
	{
		if(isWobble(ct))
		{
			vtWobbles.push_back(ct);
			vfct.push_back(makeFullClaRNATuple(ct));
			//++numberOfCanBasePairs;
		}
	}

	auto vfct2 = makeReciprocalInteractions(vfct);
	vftWobbles = finalBasePairCheck(vfct2);
	vftWobbles = addNonInteractedEdgesForWobble(vfct_formPDB, vfct2);
	checkClaRNAwithPDB(vftWobbles, vfct_formPDB);
}

void CLARNA::FindInteractions::sepNonCanBasePairs(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB)
{
	
	std::vector<FullClaRNATuple> vfct;
	// iterate all over the vtInteractions
	for(const auto& ct : vct)
	{
		if(!isCanonical(ct) && !isStacking(ct))
		{
			vtCans.push_back(ct);
			vfct.push_back(makeFullClaRNATuple(ct));
			//++numberOfCanBasePairs;
		}
	}
	
	auto vfct2 = makeReciprocalInteractions(vfct);
	vftNonCans = finalBasePairCheck(vfct2);
	vftNonCans = addNonInteractedEdgesForNonCans(vfct_formPDB, vfct2);
	checkClaRNAwithPDB(vftNonCans, vfct_formPDB);
}

void CLARNA::FindInteractions::sepBasePairs(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB)
{
	
	std::vector<FullClaRNATuple> vfct;
	// iterate all over the vtInteractions
	for(const auto& ct : vct)
	{
		if(!isStacking(ct))
		{
			vtCans.push_back(ct);
			vfct.push_back(makeFullClaRNATuple(ct));
			//++numberOfCanBasePairs;
		}
	}
	
	auto vfct2 = makeReciprocalInteractions(vfct);
	vftBasePairs= finalBasePairCheck(vfct2);
	vftBasePairs = addNonInteractedEdgesForBasePairs(vfct_formPDB, vfct2);
	checkClaRNAwithPDB(vftBasePairs, vfct_formPDB);
}

void CLARNA::FindInteractions::sepStacking(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB)
{
	std::vector<FullClaRNATuple> vfct;
	// iterate all over the vtInteractions
	for(const auto& ct : vct)
	{
		if(isStacking(ct))
		{
			vtCans.push_back(ct);
			vfct.push_back(makeFullClaRNATuple(ct));
			//++numberOfCanBasePairs;
		}
	}
	
	auto vfct2 = makeReciprocalInteractions(vfct);
	vftStacks = finalBasePairCheck(vfct2);
	vftStacks = addNonInteractedEdgesForStacks(vfct_formPDB, vfct2);
	checkClaRNAwithPDB(vftStacks, vfct_formPDB);
}


// finalBasePairCheck
std::vector<FullClaRNATuple> CLARNA::FindInteractions::finalBasePairCheck(const std::vector<FullClaRNATuple>& vfct)
{
	std::vector<FullClaRNATuple> cleanVfct;
	
	for (const auto& fct : vfct)
	{
		if(std::get<3>(fct) == ">" || std::get<3>(fct) == "<")
		{
			cleanVfct.push_back(finalStackingCheck(vfct, fct));
		}
		else
		{
			cleanVfct.push_back(finalBasePairCheck(vfct, fct));
		}
	}
	
	
	return cleanVfct;
}

FullClaRNATuple CLARNA::FindInteractions::finalBasePairCheck(const std::vector<FullClaRNATuple>& vfct, const FullClaRNATuple& fct)
{
	FullClaRNATuple cleanFct;
	auto [firstChainID, firstNuclNumber, firstNuclName, firstEdge, secondChainID, secondNuclNumber, secondNuclName, secondEdge, type, weight] = fct;
	
	// using find algorithm for see if there is corresponding reciprocal interactions

	auto itReciprocalBasePair = std::find_if(cbegin(vfct), cend(vfct), [&fct](FullClaRNATuple fctLambda)
																	{
																		return(std::get<0>(fct) == std::get<4>(fctLambda) &&
																			   std::get<1>(fct) == std::get<5>(fctLambda) &&
																			   std::get<2>(fct) == std::get<6>(fctLambda) &&
																			   std::get<3>(fct) == std::get<7>(fctLambda) &&
																			   std::get<4>(fct) == std::get<0>(fctLambda) &&
																			   std::get<5>(fct) == std::get<1>(fctLambda) &&
																			   std::get<6>(fct) == std::get<2>(fctLambda) &&
																			   std::get<7>(fct) == std::get<3>(fctLambda) &&
																			   std::get<8>(fct) == std::get<8>(fctLambda));
																	});
																	
	if(itReciprocalBasePair != cend(vfct))
	{
		cleanFct = std::make_tuple(firstChainID, firstNuclNumber, firstNuclName, firstEdge, secondChainID, secondNuclNumber, secondNuclName, secondEdge, type, weight);
	}
	else
	{
		cleanFct = std::make_tuple(firstChainID, firstNuclNumber, firstNuclName, firstEdge, "_", -1, "_", "_", "_", 0.0);
	}
	
	return cleanFct;
}

// finalStackingCheck
std::vector<FullClaRNATuple> CLARNA::FindInteractions::finalStackingCheck(const std::vector<FullClaRNATuple>& vfct)
{
	std::vector<FullClaRNATuple> cleanVfct;
	for (const auto& fct : vfct)
	{
		cleanVfct.push_back(finalStackingCheck(vfct, fct));
	}
	
	
	return cleanVfct;
}

FullClaRNATuple CLARNA::FindInteractions::finalStackingCheck(const std::vector<FullClaRNATuple>& vfct, const FullClaRNATuple& fct)
{
	FullClaRNATuple cleanFct;
	auto [firstChainID, firstNuclNumber, firstNuclName, firstEdge, secondChainID, secondNuclNumber, secondNuclName, secondEdge, type, weight] = fct;
	auto itReciprocalBasePair = std::find_if(cbegin(vfct), cend(vfct), [&fct, this](FullClaRNATuple fctLambda)
																	{
																		return(std::get<0>(fct) == std::get<4>(fctLambda) &&
																			   std::get<1>(fct) == std::get<5>(fctLambda) &&
																			   std::get<2>(fct) == std::get<6>(fctLambda) &&
																			   std::get<3>(fct) == stackingReciprocal(std::get<7>(fctLambda)) &&
																			   std::get<4>(fct) == std::get<0>(fctLambda) &&
																			   std::get<5>(fct) == std::get<1>(fctLambda) &&
																			   std::get<6>(fct) == std::get<2>(fctLambda) &&
																			   std::get<7>(fct) == stackingReciprocal(std::get<3>(fctLambda)) &&
																			   std::get<8>(fct) == std::get<8>(fctLambda));
																	});
	if(itReciprocalBasePair != cend(vfct))
	{
		cleanFct = fct;
	}
	else
	{
		cleanFct = std::make_tuple(firstChainID, firstNuclNumber, firstNuclName, firstEdge, "_", -1, "_", "_", "_", 0.0);
	}
	
	return cleanFct;
	
}


// checkClaRNAwithPDB
void CLARNA::FindInteractions::checkClaRNAwithPDB(const std::vector<FullClaRNATuple>& vfct_clarna, const std::vector<FullClaRNATuple>& vfct_pdb)
{
	// this vfct_clarna must contain non-interacted edges, and be cleaned.
	// the size of vfct_clarna and vfct_pdb must be the same
	std::string clarna_seq;
	//int num_error { 0 };
	
	if(vfct_clarna.size() != vfct_pdb.size())
	{
		cout << "error: the size of the arguments are not the same" << endl;
		exit(EXIT_FAILURE);
	}

	auto vfct_size = vfct_clarna.size();
	

	for(size_t i { 0 }; i < vfct_size; i+=m_number_involved_faces_edges)
	{
		clarna_seq += std::get<2>(vfct_clarna[i]);
		if(std::get<0>(vfct_clarna[i]) != std::get<0>(vfct_clarna[i + m_number_involved_faces_edges]) )
		{
			clarna_seq += " ";
		}
	
	}	
	
	if(clarna_seq != seqWithSeparateChains)
	{
		cout << "Choosing wrong inputs, check if PDB file and ClaRNA output file are compatible\n";
		cout << "ClaRNA: " <<  clarna_seq << endl;
		cout << "PDB:    " << seqWithSeparateChains << endl;
	}
}

/*ClaRNAMatrix CLARNA::FindInteractions::v_ClaRNATuple2ClaRNAMatrix(const std::vector<ClaRNATuple>& vct) const
{
	ClaRNAMatrix cMat;
	initClaRNAMatrix(cMat, m_sequence.size());
	for(auto const& ct : vct)
	{   
		//       0                  1               2              3                4               5                 6             7
		//auto[firstChainID, firstNuclNumber, firstNuclName, secondChainID, secondNuclNumber, secondNuclName, typeOfInteraction, weight] = ct;
		auto typeOfInteraction = std::get
		auto firstChainID = std::get<0>(fct);
		auto secondChainID = std::get<3>(fct);
		auto firstNuclNumber = std::get<1>(fct);
		auto secondNuclNumber = std::get<4>(fct);
		
		auto itID1 = std::find(possibleChainID.begin(), possibleChainID.end(), firstChainID[0]);
		auto firstChainID_num = itID1 - possibleChainID.begin();
		auto itID2 = std::find(possibleChainID.begin(), possibleChainID.end(), secondChainID[0]);
		auto secondChainID_num = itID2 - possibleChainID.begin();
		
		// W edge
		if(typeOfInteraction == "WW_cis")
		{
			std::get<0>(std::get<0>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber])) = 1;
		}
		
		if(firstEdge == "W" && secondEdge == "W" && type == "t")
		{
			std::get<0>(std::get<0>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber])) = 1;
		}
		
		if(firstEdge == "W" && secondEdge == "H" && type == "c")
		{
			std::get<0>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "W" && secondEdge == "H" && type == "t")
		{
			std::get<0>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "W" && secondEdge == "S" && type == "c")
		{
			std::get<0>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "W" && secondEdge == "S" && type == "t")
		{
			std::get<0>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		// H edge
		if(firstEdge == "H" && secondEdge == "W" && type == "c")
		{
			std::get<1>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "H" && secondEdge == "W" && type == "t")
		{
			std::get<1>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "H" && secondEdge == "H" && type == "c")
		{
			std::get<1>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "H" && secondEdge == "H" && type == "t")
		{
			std::get<1>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "H" && secondEdge == "S" && type == "c")
		{
			std::get<1>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "H" && secondEdge == "S" && type == "t")
		{
			std::get<1>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		
		// S edge
		if(firstEdge == "S" && secondEdge == "W" && type == "c")
		{
			std::get<2>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "S" && secondEdge == "W" && type == "t")
		{
			std::get<2>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "S" && secondEdge == "H" && type == "c")
		{
			std::get<2>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "S" && secondEdge == "H" && type == "t")
		{
			std::get<2>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "S" && secondEdge == "S" && type == "c")
		{
			std::get<2>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "S" && secondEdge == "S" && type == "t")
		{
			std::get<2>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
				
		// > edge
		if(firstEdge == ">" && secondEdge == ">")
		{
			std::get<3>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == ">" && secondEdge == "<")
		{
			std::get<3>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
				
		// < edge
		if(firstEdge == "<" && secondEdge == ">")
		{
			std::get<4>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
		
		if(firstEdge == "<" && secondEdge == "<")
		{
			std::get<4>(cMat[firstChainID_num * firstNuclNumber][secondChainID_num * secondNuclNumber]) = 1;
		}
	}
	
	
	for(const auto& vec : cMat)
	{
		for(const auto& elem : vec)
			stdExt::tuplePrint2(elem);
		cout << endl;
	}
	return cMat;
}*/
