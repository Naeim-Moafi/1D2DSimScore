#include "FindInteractions.h"
//#include "../format.h"
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

using std::cout;
using std::endl;
using std::cerr;




// constructor
// default constructor
SS::FindInteractions::FindInteractions()
{
	m_numberOfCanBasePairs = 0;
	m_numberOfNonCanBasePairs = 0;
	m_numberOfAllBasePairs = 0;
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

// separateOpensDotsCloses
OpenDotCloseIndices SS::FindInteractions::separateOpensDotsCloses(BasePairSymbols basePairSymbols, std::string ss) const
{
	std::vector<int> vOpens, vDots, vCloses;
	// for dot we only need to use ss( secondary structure inside input files)
	for(size_t i{0}; i < ss.size(); ++i)
	{
		if(ss[i] == '.')
		{
			// the indices in real world start form one :)
			vDots.push_back(i + 1);
		}
	}
	
	// BasePairSymbols maps open and close indices we check for these pairs inside the ss 
	// for each character one by one
	for(size_t i{0}; i < ss.size(); ++i)
	{
		auto itOpen = std::find_if(std::begin(basePairSymbols), std::end(basePairSymbols), [&ss, i](BasePairSymbol bpChar){return (bpChar.first == ss[i]);});
		auto itClose = std::find_if(std::begin(basePairSymbols), std::end(basePairSymbols), [&ss, i](BasePairSymbol bpChar){return (bpChar.second == ss[i]);});
		
		if(itOpen != std::end(basePairSymbols))
		{
			vOpens.push_back(i + 1);
		}
		
		if(itClose != std::end(basePairSymbols))
		{
			vCloses.push_back( i + 1);
		}
	}
	
	if (vOpens.size() != vCloses.size())
	{
		throw std::invalid_argument("The number of open and close are not equal");
	}
	
	return OpenDotCloseIndices{ .vOpens = vOpens, .vDots = vDots, .vCloses = vCloses};
}

// ss2IndicesPair
SSMap SS::FindInteractions::ss2IndicesPair(BasePairSymbols basePairSymbols, std::string ss)
{
	SSMap basePairsMap;
	OpenDotCloseIndices odcIndices;
	try
	{
		odcIndices = separateOpensDotsCloses(basePairSymbols, ss);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	// insert non interacted nucleotides
	for(const auto& dotIndex : odcIndices.vDots)
	{
		basePairsMap.insert(std::make_pair(dotIndex, -1));
	}
	
	int numberOfExpectedInteractions = odcIndices.vOpens.size();;
	int numberOfDetectedInteractions = 0;
	
	// finding the base-pairs with iterating over the open and close symbols vectors
	for(int  i = odcIndices.vOpens.size() - 1; i >= 0; --i)
	{
		for(size_t j{0}; j < odcIndices.vCloses.size(); ++j)
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
				basePairsMap.insert(std::make_pair(odcIndices.vOpens[i], odcIndices.vCloses[j]));
				basePairsMap.insert(std::make_pair(odcIndices.vCloses[j], odcIndices.vOpens[i]));
				++m_numberOfAllBasePairs;
				ss.replace(ss.begin() + odcIndices.vOpens[i] - 1, ss.begin() + odcIndices.vOpens[i], ".");
				ss.replace(ss.begin() + odcIndices.vCloses[j] - 1, ss.begin() + odcIndices.vCloses[j], ".");
				++numberOfDetectedInteractions;
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
	
	return basePairsMap;
}

// extractInteraction
SSMap SS::FindInteractions::extractInteractions(std::string ss)
{
	BasePairSymbols basePairSymbols;
	SSMap ssMap;
	
	// feed the basePairSymbols with information 
	// or in other words, make dictionary of open and close Symbols
	for(size_t i{0}; i < openSymbols.size(); ++i)
	{
		basePairSymbols.insert(std::make_pair(openSymbols[i], closeSymbols[i]));
	}
	
	try
	{
		ssMap = ss2IndicesPair( basePairSymbols, ss);
	}
	catch (const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	return ssMap; 
}

// readInputFile
SSMap SS::FindInteractions::readInputFile(std::filesystem::path path)
{
	SSMap ssMap;
	//cout << fmt::format("Start reading file {}", path.filename().string()) << endl;
	//cout << "Start reading file " << path.filename().string() << endl;
	
	int numberOfLine = 0;
	std::ifstream inpFile(path.string().c_str());
	std::string ss;
	
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
		std::getline(inpFile, ss);
		++numberOfLine;
		if(ss.size() == 0)
		{ // break for emty line
			break;
		}
		
		if(numberOfLine > 1)
		{
			throw std::invalid_argument("the format of ss is not correct");
		}
		
		try{ 
			//cout << "The secondary structure: " <<  ss << endl;
			ssSize = ss.size();
			ssMap = extractInteractions(ss);
		}
		catch(const std::invalid_argument& ex)
		{
			cerr << ex.what() << endl;
			//cout << fmt::format("Please check file '{}' in follwing directory\n'{}'", path.filename().string(), path.parent_path().string()) << endl;
			cout << "Please check the file '"<< path.filename().string() <<"' in the followng directory\n";
			cout << "'" << path.relative_path().string() << "'" << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	inpFile.close();
	
	return ssMap;
}

// readSeqFile
std::vector<std::string> SS::FindInteractions::readSeqFile(std::filesystem::path path)
{
	// sequence of each chain will be save as differne elements of the following vector
	std::vector<std::string> vSeqOfChains;
	std::ifstream inpFile(path.string().c_str());
	std::string seq;
	std::ostringstream ossSeq;
	
	
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
	
	
	// defaule separator of the operator >> is space
	// in this way we can find the seqence of different chains
	while(inpFile >> seq)
	{
		vSeqOfChains.push_back(seq);
	}
	
	// copy whol sequence in a string
	copy(std::begin(vSeqOfChains), std::end(vSeqOfChains), std::ostream_iterator<std::string>(ossSeq, " "));
	seqWithSeparateChains = ossSeq.str();
	
	inpFile.close();
	return vSeqOfChains;	
}

// isCanonical
bool SS::FindInteractions::isCanonical(std::vector<std::string> vSeqOfChains, std::pair<int, int> ssPair) const
{	
	if(ssPair.second != -1)
	{
		if((vSeqOfChains[0][ssPair.first - 1] == 'A' && vSeqOfChains[0][ssPair.second - 1] == 'U') 
		   || (vSeqOfChains[0][ssPair.first - 1] == 'U' && vSeqOfChains[0][ssPair.second - 1] == 'A') )
		{
			return true;
		}
		
		if((vSeqOfChains[0][ssPair.first - 1] == 'C' && vSeqOfChains[0][ssPair.second - 1] == 'G')
		   || (vSeqOfChains[0][ssPair.first - 1] == 'G' && vSeqOfChains[0][ssPair.second - 1] == 'C'))
		{
			return true;
		}
		
		if(m_isWobble_canonical)
		{
			if((vSeqOfChains[0][ssPair.first - 1] == 'G' && vSeqOfChains[0][ssPair.second - 1] == 'U')
			   || (vSeqOfChains[0][ssPair.first - 1] == 'U' && vSeqOfChains[0][ssPair.second - 1] == 'G'))
			{
				return true;
			}
		}
	}
	
	return false;
}

bool SS::FindInteractions::isCanonical(std::string seq, std::pair<int, int> ssPair) const
{	
	if(ssPair.second != -1)
	{
		if((seq[ssPair.first - 1] == 'A' && seq[ssPair.second - 1] == 'U') 
		   || (seq[ssPair.first - 1] == 'U' && seq[ssPair.second - 1] == 'A') )
		{
			return true;
		}
		
		if((seq[ssPair.first - 1] == 'C' && seq[ssPair.second - 1] == 'G')
		   || (seq[ssPair.first - 1] == 'G' && seq[ssPair.second - 1] == 'C'))
		{
			return true;
		}
		
		if(m_isWobble_canonical)
		{
			if((seq[ssPair.first - 1] == 'G' && seq[ssPair.second - 1] == 'U')
			   || (seq[ssPair.first - 1] == 'U' && seq[ssPair.second - 1] == 'G'))
			{
				return true;
			}
		}
		
	}
	return false;
}

// mapSS2seq
void SS::FindInteractions::mapSS2seq(std::vector<std::string> vSeqOfChains, SSMap ssMap)
{
	
	if(vSeqOfChains.size() > 1)
	{
		//cout << fmt::format("{} feature is not ready for the vector size {} yet", __func__, vSeqOfChains.size()) << endl;
		cout << __func__ <<  " feature is not ready for the vector size" <<  vSeqOfChains.size() << "yet\n";
		exit(EXIT_FAILURE);
	}
	
	for(const auto& ssPair : ssMap)
	{
		if(isCanonical(vSeqOfChains, ssPair))
		{
			ssMap_can.insert(std::make_pair(ssPair.first, ssPair.second));
			ssMap_can.insert(std::make_pair(ssPair.second, ssPair.first));
			ssMap_noncan.insert(std::make_pair(ssPair.first, -1));
			++m_numberOfCanBasePairs;
		}
		else
		{
			ssMap_can.insert(std::make_pair(ssPair.first, -1));
			ssMap_noncan.insert(std::make_pair(ssPair.first, ssPair.second));
			if(ssPair.second != -1)
			{
				
				ssMap_noncan.insert(std::make_pair(ssPair.second, ssPair.first));
				//cout << ssPair.first << " " << ssPair.second << endl;
				++m_numberOfNonCanBasePairs;
			}
		}
	}
}


void SS::FindInteractions::mapSS2seq(std::string seq, SSMap ssMap)
{	
	for(const auto& ssPair : ssMap)
	{
		if(isCanonical(seq, ssPair))
		{
			ssMap_can.insert(std::make_pair(ssPair.first, ssPair.second));
			ssMap_can.insert(std::make_pair(ssPair.second, ssPair.first));
			ssMap_noncan.insert(std::make_pair(ssPair.first, -1));
			++m_numberOfCanBasePairs;
		}
		else
		{
			ssMap_can.insert(std::make_pair(ssPair.first, -1));
			ssMap_noncan.insert(std::make_pair(ssPair.first, ssPair.second));
			if(ssPair.second != -1)
			{
				ssMap_noncan.insert(std::make_pair(ssPair.second, ssPair.first));
				//cout << ssPair.first << " " << ssPair.second << endl;
				++m_numberOfNonCanBasePairs;
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
	SS::FindInteractions fi2(true, true);
	cout << fmt::format("fi: isWobble_canonical --> {}, withSeq- --> {}", fi.get_isWobble_canonical(), fi.get_withSeq()) << endl;
	cout << fmt::format("fi2: isWobble_canonical --> {}, withSeq- --> {}", fi2.get_isWobble_canonical(), fi2.get_withSeq()) << endl;
	//fi.set_isWobble_canonical(true);
	fi.set_withSeq(true);
	fi2.set_withSeq(false);
	cout << "after set" << endl;
	cout << fmt::format("fi: isWobble_canonical --> {}, withSeq- --> {}", fi.get_isWobble_canonical(), fi.get_withSeq()) << endl;
	cout << fmt::format("fi2: isWobble_canonical --> {}, withSeq- --> {}", fi2.get_isWobble_canonical(), fi2.get_withSeq()) << endl;
	
	
	SSMap ssMap = fi.readInputFile("dotBracketRef.SS");
	
	for(const auto& ssPair : ssMap)
	{
		cout << fmt::format("({:2}, {:2})", ssPair.first, ssPair.second) << endl;
	}
	
	vector<string> vSeqOfChains = fi.readSeqFile("SeqForDotBracket.seq");	
	//vector<string> vSeqOfChains = fi.readSeqFile("Seq.seq");
	cout << fi.seqWithSeparateChains << endl;
	
	if(fi.get_withSeq())
	{
		fi.mapSS2seq(vSeqOfChains, ssMap);
		cout << fmt::format("#Canonical interaction: {}", fi.get_numberOfCanBasePairs()) << endl;
		for(const auto& ssPair : fi.ssMap_can)
		{
			cout << fmt::format("({:2}, {:2})", ssPair.first, ssPair.second) << endl;
		}
		
		cout << endl;
		cout << fmt::format("#Non-canonical interaction: {}", fi.get_numberOfNonCanBasePairs()) << endl;
		for(const auto& ssPair : fi.ssMap_noncan)
		{
			cout << fmt::format("({:2}, {:2})", ssPair.first, ssPair.second) << endl;
		}
	}
	
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
