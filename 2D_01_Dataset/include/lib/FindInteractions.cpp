#include "FindInteractions.h"
//#include "../format.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <chrono>

//using namespace std;
using namespace std::chrono;

using std::cout;
using std::endl;
using std::cerr;

void initSSMatrix(SSMatrix& ssMat, size_t length)
{
	ssMat.clear();
	for(size_t i { 0 }; i < length; ++i)
	{
		std::vector<int> v_tmp;
		for(size_t j { 0 }; j < length; ++j)
		{
			v_tmp.push_back(0);
		}
		ssMat.emplace_back(v_tmp);
	}
}


// constructor
// default constructor
SS::FindInteractions::FindInteractions()
{
	m_numberOfCanBasePairs = 0;
	m_numberOfNonCanBasePairs = 0;
	m_numberOfAllBasePairs = 0;
	m_isWobble_canonical = false;
	m_withSeq = false;
	ssSize = 0;
}

SS::FindInteractions::FindInteractions(bool isWobble_canonical, bool withSeq)
	: SS::FindInteractions()
{
	m_isWobble_canonical = isWobble_canonical;
	m_withSeq = withSeq;
}

// setter and getter of the private data members
int SS::FindInteractions::get_numberOfCanBasePairs() const
{
	return m_numberOfCanBasePairs/2;
}

int SS::FindInteractions::get_numberOfNonCanBasePairs() const
{
	return m_numberOfNonCanBasePairs/2;
}

int SS::FindInteractions::get_numberOfAllBasePairs() const
{
	return m_numberOfAllBasePairs/2;
}

bool SS::FindInteractions::get_isWobble_canonical() const
{
	return m_isWobble_canonical;
}

void SS::FindInteractions::set_isWobble_canonical(bool isWobble_canonical)
{
	m_isWobble_canonical = isWobble_canonical;
}

bool SS::FindInteractions::get_withSeq() const
{
	return m_withSeq;
}

void SS::FindInteractions::set_withSeq(bool withSeq)
{
	m_withSeq = withSeq;
}

void SS::FindInteractions::init_sequence(const std::filesystem::path& seqPath)
{
	std::ifstream seqFile(seqPath.string());
	std::string chainSeq;
	std::vector<std::string> v_chainsSeq;
	while(seqFile >> chainSeq)
	{
		v_chainsSeq.emplace_back(chainSeq);
		m_sequence +=  chainSeq;
	}

	// separating the sequence with its chains;
	std::ostringstream oss;
	std::copy(v_chainsSeq.begin(), v_chainsSeq.end(), std::ostream_iterator<std::string>(oss, ""));
	m_seqWithSeparateChains = oss.str();
}

// separateOpensDotsCloses
OpenDotCloseIndices SS::FindInteractions::separateOpensDotsCloses(const BasePairSymbols& basePairSymbols, const std::vector<std::string>& v_ss) const
{
	std::vector<int> vOpens, vDots, vCloses, vBrackets;
	
	for( const auto& ss : v_ss)
	{
		for(size_t i { 0 }; i < ss.size(); ++i)
		{
			auto itOpen = std::find_if(std::begin(basePairSymbols), std::end(basePairSymbols), [&ss, i](BasePairSymbol bpChar){return (bpChar.first == ss[i]);});
			auto itClose = std::find_if(std::begin(basePairSymbols), std::end(basePairSymbols), [&ss, i](BasePairSymbol bpChar){return (bpChar.second == ss[i]);});
			
			if(itOpen != std::end(basePairSymbols))
			{
				vOpens.push_back(i + 1);
				vBrackets.push_back(i + 1);
			}
			
			if(itClose != std::end(basePairSymbols))
			{
				vCloses.push_back(i + 1);
				vBrackets.push_back(i + 1);
			}
		}
	}
	// finding the dots
	for(size_t i { 0 }; i < v_ss[0].size(); ++i)
	{
		auto itBracket = find(vBrackets.begin(), vBrackets.end(), i + 1);
		if(itBracket == vBrackets.end())
		{
			vDots.push_back(i+1);
		}
	}
	
	if (vOpens.size() != vCloses.size())
	{
		throw std::invalid_argument("The number of open and close are not equal");
	}
	
	return OpenDotCloseIndices{.vOpens = vOpens, .vDots = vDots, .vCloses = vCloses};
}

// ss2IndicesPair
SSContacts SS::FindInteractions::ss2IndicesPair(BasePairSymbols& basePairSymbols, std::vector<std::string> v_ss)
{
	SSContacts ssc;
	OpenDotCloseIndices odcIndices;
	try
	{
		odcIndices = separateOpensDotsCloses(basePairSymbols, v_ss);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	// insert non interacted nucleotides
	for(const auto& dotIndex : odcIndices.vDots)
	{
		ssc.emplace_back(std::make_pair(dotIndex, -1));
	}
	
	int numberOfExpectedInteractions = odcIndices.vOpens.size();;
	int numberOfDetectedInteractions = 0;
	for(auto& ss : v_ss)
	{
		// finding the base-pairs with iterating over the open and close symbols vectors
		for(int  i = odcIndices.vOpens.size() - 1; i >= 0; --i)
		{
			for(size_t j { 0 }; j < odcIndices.vCloses.size(); ++j)
			{
				// for base paring the index of the close symbol always must be greater than open one
				if(odcIndices.vCloses[j] < odcIndices.vOpens[i])
				{
					continue; // do not consider it
				}
				
				// basePairSymbols is a map which has keyType = char, and valueType = char
				// so we get key from ss input file and compare the value (close symbol) to the one in the ss input file.
				if(basePairSymbols[ss[odcIndices.vOpens[i] - 1]] == ss[odcIndices.vCloses[j] - 1])
				{
					ssc.emplace_back(std::make_pair(odcIndices.vOpens[i], odcIndices.vCloses[j]));
					ssc.emplace_back(std::make_pair(odcIndices.vCloses[j], odcIndices.vOpens[i]));
					++m_numberOfAllBasePairs;
					ss.replace(ss.begin() + odcIndices.vOpens[i] - 1, ss.begin() + odcIndices.vOpens[i], ".");
					ss.replace(ss.begin() + odcIndices.vCloses[j] - 1, ss.begin() + odcIndices.vCloses[j], ".");
					++numberOfDetectedInteractions;
				}
			}
		}
	}
	// checking the number of interactions
	if(numberOfExpectedInteractions != numberOfDetectedInteractions)
	{
		//cout << fmt::format("the number of expected base pairs is {} but {} number of base pairs are found", numberOfExpectedInteractions, numberOfDetectedInteractions) << endl;
		cout << "the number of expected base pairs is" << numberOfExpectedInteractions << "but" <<  numberOfDetectedInteractions << "number of base pairs are found" << endl;
		throw std::invalid_argument("Something wrong in the secondary structure");
	}
	
	
	return ssc;
}

SSMap SS::FindInteractions::SSContacts2SSMap(const SSContacts& ssc)
{
	SSMap ssMap;
	for(const auto& pair : ssc)
	{
		ssMap.insert(pair);
	}
	
	return ssMap;
}

SSMatrix SS::FindInteractions::SSContacts2SSMatrix(const SSContacts& ssc)
{
	SSMatrix ssmat;
	if(m_sequence.size() != 0)
	{
		initSSMatrix(ssmat, m_sequence.size());
	}
	else
	{
		initSSMatrix(ssmat, m_SSMap.size());
	}
	for(const auto& pair : ssc)
	{
		if(pair.second != -1)
		{
			ssmat[pair.first - 1][pair.second - 1] = 1;
		}
	}
	return ssmat;
}

SSMatrix SS::FindInteractions::extractInteractions_matrix(const std::vector<std::string>& v_ss)
{
	SSMatrix ssMatrix;
	SSContacts ssContacts;
	// make dictionary of the open and correspoding bracket
	BasePairSymbols basePairSymbols;
	for(size_t i{0}; i < openSymbols.size(); ++i)
	{
		basePairSymbols.insert(std::make_pair(openSymbols[i], closeSymbols[i]));
	}
	try
	{
		ssContacts = ss2IndicesPair(basePairSymbols, v_ss);
		ssMatrix = SSContacts2SSMatrix(ssContacts);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	return ssMatrix;
}

SSMap SS::FindInteractions::extractInteractions_vector(const std::vector<std::string>& v_ss)
{
	SSMap ssMap;
	SSContacts ssContacts;
	// make dictionary of the open and correspoding bracket
	BasePairSymbols basePairSymbols;
	for(size_t i{0}; i < openSymbols.size(); ++i)
	{
		basePairSymbols.insert(std::make_pair(openSymbols[i], closeSymbols[i]));
	}
	
	try
	{
		ssContacts = ss2IndicesPair(basePairSymbols, v_ss);
		ssMap = SSContacts2SSMap(ssContacts);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	return ssMap;
}

// init_structure
void SS::FindInteractions::init_structure(const std::filesystem::path& path)
{
	//cout << fmt::format("Start reading file {}", path.filename().string()) << endl;
	//cout << "Start reading file " << path.filename().string() << endl;
	
	//int numberOfLine = 0;
	std::ifstream inpFile(path.string().c_str());
	std::string line;
	std::vector<std::string> v_ss;
	
	if(inpFile.fail())
	{
		if(!path.parent_path().empty())
		{
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), path.relative_path().string()) << endl;
			cout << "Please check the file '"<< path.filename().string() <<"' in the followng directory\n";
			cout << "'" << path.relative_path().string() << "'" << endl;
		}
		else
		{
			std::filesystem::path p = std::filesystem::current_path();
			//cout << fmt::format("Please check the file '{}' in the followng directory\n'{}'", path.filename().string(), p.relative_path().string()) << endl;
			cout << "Please check the file '"<< path.filename().string() <<"' in the followng directory\n";
			cout << "'" << p.relative_path().string() << "'" << endl;
		}
		exit(EXIT_FAILURE);
	}

	while(!inpFile.fail())
	{
		std::getline(inpFile, line);
		if(line.size() == 0)
		{ // break for emty line
			break;
		}
		
		std::istringstream iss(line);
		std::string ss;
		std::string chainSS;
		while(iss >> chainSS)
		{
			ss += chainSS;
		}
		
		v_ss.emplace_back(ss);
		ssSize = ss.size();
	}
	
		//std::copy(v_chainsSeq.begin(), v_chainsSeq.end(), std::ostream_iterator<std::string>(cout, " "));
	
	try{ 
		//cout << "The secondary structure: " <<  ss << endl;
		m_SSMap = extractInteractions_vector(v_ss);
		m_Matrix_all = extractInteractions_matrix(v_ss);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		//cout << fmt::format("Please check file '{}' in follwing directory\n'{}'", path.filename().string(), path.parent_path().string()) << endl;
		cout << "Please check the file '"<< path.filename().string() <<"' in the followng directory\n";
		cout << "'" << path.relative_path().string() << "'" << endl;
		exit(EXIT_FAILURE);
	}
	
	inpFile.close();
}


// isCanonical
bool SS::FindInteractions::isCanonical(std::pair<int, int> ssPair) const
{	
	if(ssPair.second != -1)
	{
		if((m_sequence[ssPair.first - 1] == 'A' && (m_sequence[ssPair.second - 1] == 'U' || m_sequence[ssPair.second - 1] == 'T')) 
		   || ((m_sequence[ssPair.first - 1] == 'U' || m_sequence[ssPair.first - 1] == 'T') && m_sequence[ssPair.second - 1] == 'A') )
		{
			return true;
		}
		
		if((m_sequence[ssPair.first - 1] == 'C' && m_sequence[ssPair.second - 1] == 'G')
		   || (m_sequence[ssPair.first - 1] == 'G' && m_sequence[ssPair.second - 1] == 'C'))
		{
			return true;
		}
		
		if(m_isWobble_canonical)
		{
			if((m_sequence[ssPair.first - 1] == 'G' && (m_sequence[ssPair.second - 1] == 'U' || m_sequence[ssPair.second - 1] == 'T'))
			   || ((m_sequence[ssPair.first - 1] == 'U' || m_sequence[ssPair.first - 1] == 'T') && m_sequence[ssPair.second - 1] == 'G'))
			{
				return true;
			}
		}
	}
	
	return false;
}

bool SS::FindInteractions::isWobble(std::pair<int, int> ssPair) const
{
	if((m_sequence[ssPair.first - 1] == 'G' && (m_sequence[ssPair.second - 1] == 'U' || m_sequence[ssPair.second - 1] == 'T'))
	   || ((m_sequence[ssPair.first - 1] == 'U' || m_sequence[ssPair.first - 1] == 'T') && m_sequence[ssPair.second - 1] == 'G'))
	{
		return true;
	}
	
	return false;
}

// mapSS2seq
void SS::FindInteractions::mapSS2seq(SSMap ssMap)
{
	initSSMatrix(m_Matrix_can, m_sequence.size());
	initSSMatrix(m_Matrix_noncan, m_sequence.size());
	initSSMatrix(m_Matrix_wobble, m_sequence.size());
	for(const auto& ssPair : ssMap)
	{
		if(isCanonical(ssPair))
		{
			m_SSMap_can.insert(std::make_pair(ssPair.first, ssPair.second));
			m_SSMap_can.insert(std::make_pair(ssPair.second, ssPair.first));			
			m_SSMap_noncan.insert(std::make_pair(ssPair.first, -1));
			if(isWobble(ssPair))
			{
				m_SSMap_wobble.insert(std::make_pair(ssPair.first, ssPair.second));
				m_SSMap_wobble.insert(std::make_pair(ssPair.second, ssPair.first));
				
				m_Matrix_wobble[ssPair.first - 1][ssPair.second - 1] = 1;
				m_Matrix_wobble[ssPair.second - 1][ssPair.first - 1] = 1;
			}
			else
			{
				m_SSMap_wobble.insert(std::make_pair(ssPair.first, -1));
			}
			// add the information into the canonical interactins' matrix
			m_Matrix_can[ssPair.first - 1][ssPair.second - 1] = 1;
			m_Matrix_can[ssPair.second - 1][ssPair.first - 1] = 1;
			++m_numberOfCanBasePairs;
		}
		else
		{
			m_SSMap_can.insert(std::make_pair(ssPair.first, -1));
			m_SSMap_noncan.insert(std::make_pair(ssPair.first, ssPair.second));
			if(ssPair.second != -1)
			{
				m_SSMap_noncan.insert(std::make_pair(ssPair.second, ssPair.first));
				if(isWobble(ssPair))
				{
					m_SSMap_wobble.insert(std::make_pair(ssPair.first, ssPair.second));
					m_SSMap_wobble.insert(std::make_pair(ssPair.second, ssPair.first));
						
					m_Matrix_wobble[ssPair.first - 1][ssPair.second - 1] = 1;
					m_Matrix_wobble[ssPair.second - 1][ssPair.first - 1] = 1;
				}
				else
				{
					m_SSMap_wobble.insert(std::make_pair(ssPair.first, -1));
				}
				
				// add the information into the canonical interactins' matrix
				m_Matrix_noncan[ssPair.first - 1][ssPair.second - 1] = 1;
				m_Matrix_noncan[ssPair.second - 1][ssPair.first - 1] = 1;
				++m_numberOfNonCanBasePairs;
			}
			else
			{
				m_SSMap_wobble.insert(std::make_pair(ssPair.first, -1));				
			}
		}
	}
}

//#define __MAIN__
#ifdef __MAIN__
using namespace std;
//using namespace SS;
using namespace filesystem;
 
int main()
{
	SS::FindInteractions fi;
	fi.init_sequence("samples/Cruciform.fasta");
	fi.init_structure("samples/Cruciform.SS");
	fi.mapSS2seq(fi.m_SSMap);
	
	cout << fi.m_Matrix_wobble.size() << endl;
	
	for(auto const& entry : fi.m_Matrix_wobble)
	{
		for(size_t i = 0; i < entry.size(); ++i)
		{
			if(entry[i] == 0)
				cout << ". ";
			else
				cout << "* ";
		}
		cout << endl;
	}
	//SS::FindInteractions fi2(true, true);
	//cout << fmt::format("fi: isWobble_canonical --> {}, withSeq- --> {}", fi.get_isWobble_canonical(), fi.get_withSeq()) << endl;
	//cout << fmt::format("fi2: isWobble_canonical --> {}, withSeq- --> {}", fi2.get_isWobble_canonical(), fi2.get_withSeq()) << endl;
	////fi.set_isWobble_canonical(true);
	//fi.set_withSeq(true);
	//fi2.set_withSeq(false);
	//cout << "after set" << endl;
	//cout << fmt::format("fi: isWobble_canonical --> {}, withSeq- --> {}", fi.get_isWobble_canonical(), fi.get_withSeq()) << endl;
	//cout << fmt::format("fi2: isWobble_canonical --> {}, withSeq- --> {}", fi2.get_isWobble_canonical(), fi2.get_withSeq()) << endl;
	
	
	//SSMap ssMap = fi.readInputFile("dotBracketRef.SS");
	
	//for(const auto& ssPair : ssMap)
	//{
		//cout << fmt::format("({:2}, {:2})", ssPair.first, ssPair.second) << endl;
	//}
	
	//vector<string> vSeqOfChains = fi.readSeqFile("SeqForDotBracket.seq");	
	////vector<string> vSeqOfChains = fi.readSeqFile("Seq.seq");
	//cout << fi.seqWithSeparateChains << endl;
	
	//if(fi.get_withSeq())
	//{
		//fi.mapSS2seq(vSeqOfChains, ssMap);
		//cout << fmt::format("#Canonical interaction: {}", fi.get_numberOfCanBasePairs()) << endl;
		//for(const auto& ssPair : fi.ssMap_can)
		//{
			//cout << fmt::format("({:2}, {:2})", ssPair.first, ssPair.second) << endl;
		//}
		
		//cout << endl;
		//cout << fmt::format("#Non-canonical interaction: {}", fi.get_numberOfNonCanBasePairs()) << endl;
		//for(const auto& ssPair : fi.ssMap_noncan)
		//{
			//cout << fmt::format("({:2}, {:2})", ssPair.first, ssPair.second) << endl;
		//}
	//}
	
	//cout <<  "testing the mismatches" << endl;
	//vector<int>  v  { 1, 2, 3, 4, 5, 6 };
	//vector<int>  v2 { 1, 7, 3, 8, 5, 0 };
	//auto vMismatches = stdExt::mismatches(begin(v), end(v), begin(v2), end(v2));
	//cout << "in main" << endl;
	//for( auto& mismatch : vMismatches)
	//{
		//cout << *mismatch.first << " " << *mismatch.second << endl;
	//}
	
	//cout << "another version" << endl;
	//auto vMismatches2 = stdExt::mismatches(begin(v), end(v), begin(v2), end(v2), [](auto e1, auto e2){return e2 > e1;});
	
	//for( auto& mismatch : vMismatches2)
	//{
		//cout << *mismatch.first << " " << *mismatch.second << endl;
	//}

}

#endif //__MAIN__
