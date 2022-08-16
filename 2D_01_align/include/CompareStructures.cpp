#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

using std::cout;
using std::endl;
using std::cerr;

const std::vector<std::string> vInteractionTypes1 = {"Cans,", "wobbles,", "NonCans,", "AllBasePairs,"};
const std::vector<std::string> vInteractionTypes2 = {"Cans + wobbles,", "wobbles,", "NonCans,", "AllBasePairs,"};
//constructors
BLAST::CompareStructures::CompareStructures()
{
	requestedInteractions = "CANS";
}

BLAST::CompareStructures::CompareStructures(std::string requestedInteractions)
: requestedInteractions(requestedInteractions) {}

// readStructures
void BLAST::CompareStructures::readFile(std::filesystem::path blastPath)
{
	findInteractions.readInputFile(blastPath);
	ref.all = findInteractions.resetIndices(findInteractions.m_vRanges[0], findInteractions.m_vSSMaps[0]);
	ref.cans = findInteractions.resetIndices(findInteractions.m_vRanges[0],  findInteractions.m_vSSMaps_can[0]);
	ref.nonCans = findInteractions.resetIndices(findInteractions.m_vRanges[0],  findInteractions.m_vSSMaps_noncan[0]);
	query.all = findInteractions.resetIndices(findInteractions.m_vRanges[1],  findInteractions.m_vSSMaps[1]);
	query.cans = findInteractions.resetIndices(findInteractions.m_vRanges[1],  findInteractions.m_vSSMaps_can[1]);
	query.nonCans = findInteractions.resetIndices(findInteractions.m_vRanges[1],  findInteractions.m_vSSMaps_noncan[1]);
	findInteractions.structure2Can();
	findInteractions.structure2NonCan();
	
}

int BLAST::CompareStructures::find_max_n_positives(const SSMap& ref_, const SSMap& query_)
{
	int np_ref = 0;
	int np_query = 0;
	for(auto pair : ref_)
	{
		if(pair.second != -1)
		{
			++np_ref;
		}
	}
	
	
	for(auto pair : query_)
	{
		if(pair.second != -1)
		{
			++np_query;
		}
	}
	
	if(np_ref > np_query)
	{
		return np_ref;
	}
	
	return np_query;
}

ConfusionMatrixTuple_blast BLAST::CompareStructures::calcConfusionMatrix()
{
	using std::cerr;
	using std::endl;
	using std::cout;
	
	int TP, TN, FP, FN, TPx, TNx, FPx, FNx;
	
	TP = TN = FP = FN = TPx = TNx = FPx = FNx = 0;
	
	cout << "All:\n";
	cout << "\tref:\n";
	findInteractions.showInformationInRange(findInteractions.m_vRanges[0], ref.all, 0);
	cout << "\tquery\n";
	findInteractions.showInformationInRange(findInteractions.m_vRanges[1], query.all, 1);
	
	
	
	for(size_t firstNucl = 1; firstNucl <=  ref.all.size(); ++firstNucl)
	{
		// the first nucleotide in the pairs are the same
		// so we compare only the second nucleotides
		// [example: vRefSS[i] gives the second nuleotide involved in each pair for the reference structure]
		if(ref.all[firstNucl] == query.all[firstNucl] && ref.all[firstNucl] != -1 )
		{
			if(ref.all[firstNucl] != -2 && ref.all[firstNucl] != -3)
			{
				++TP;
			}
			else
			{
				++TPx;
			}
		}
		
		if(ref.all[firstNucl] != query.all[firstNucl])
		{
			if(ref.all[firstNucl] > 0)
			{
				switch(query.all[firstNucl])
				{
					case -1: ++FN; break;
					case -2: [[fallthrough]];
					case -3: ++FPx;
					default: ++FP; break;
				}
			}
			
			if(ref.all[firstNucl] == -1)
			{
				switch(query.all[firstNucl])
				{
					case -2 : ++FPx; break;
					case -3 : ++FPx; break;
					default : ++FP; break;
				}
			}
			
			if(ref.all[firstNucl] == -2)
			{
				switch(query.all[firstNucl])
				{
					case -1 : ++FNx; break;
					default : ++FPx; break;
				}
			}
			
			if(ref.all[firstNucl] == -3)
			{
				switch(query.all[firstNucl])
				{
					case -1 : ++FNx; break;
					default : ++FPx; break;
				}
			}
		}
		
		if(ref.all[firstNucl] == query.all[firstNucl] && ref.all[firstNucl] == -1)
		{
			++TN;
		}
		
	}
	
	return std::make_tuple(TP, TN, FP, FN, TPx, TNx, FPx, FNx);
	
}

ConfusionMatrixTuple_blast BLAST::CompareStructures::calcConfusionMatrixCans()
{
	using std::cerr;
	using std::endl;
	using std::cout;
	
	int TP, TN, FP, FN, TPx, TNx, FPx, FNx;
	
	TP = TN = FP = FN = TPx = TNx = FPx = FNx = 0;
	auto itE = std::find_if(requestedInteractions.begin(), requestedInteractions.end(), [](char ch){return (ch == 'E' || ch == 'e');});
	if(itE != requestedInteractions.end())
	{
		cout << "Cans+Wobble:\n";
		findInteractions.set_isWobble_canonical(true);
	}
	else
	{
		cout << "Cans:\n";
	}
	
	cout << "\tref:\n";
	findInteractions.showInformationInRange(findInteractions.m_vRanges[0], ref.cans, 0, "c");
	cout << "\tquery\n";
	findInteractions.showInformationInRange(findInteractions.m_vRanges[1], query.cans, 1, "c");
	
	for(size_t firstNucl = 1; firstNucl <=  ref.all.size(); ++firstNucl)
	{
		// the first nucleotide in the pairs are the same
		// so we compare only the second nucleotides
		// [example: vRefSS[i] gives the second nuleotide involved in each pair for the reference structure]
		if(ref.cans[firstNucl] == query.cans[firstNucl] && ref.cans[firstNucl] != -1 )
		{
			if(ref.cans[firstNucl] != -2 && ref.cans[firstNucl] != -3)
			{
				++TP;
			}
			else
			{
				++TPx;
			}
		}
		
		if(ref.cans[firstNucl] != query.cans[firstNucl])
		{
			if(ref.cans[firstNucl] > 0)
			{
				switch(query.cans[firstNucl])
				{
					case -1: ++FN; break;
					case -2: [[fallthrough]];
					case -3: ++FPx;
					default: ++FP; break;
				}
			}
			
			if(ref.cans[firstNucl] == -1)
			{
				switch(query.cans[firstNucl])
				{
					case -2 : ++FPx; break;
					case -3 : ++FPx; break;
					default : ++FP; break;
				}
			}
			
			if(ref.cans[firstNucl] == -2)
			{
				switch(query.cans[firstNucl])
				{
					case -1 : ++FNx; break;
					default : ++FPx; break;
				}
			}
			
			if(ref.cans[firstNucl] == -3)
			{
				switch(query.cans[firstNucl])
				{
					case -1 : ++FNx; break;
					default : ++FPx; break;
				}
			}
		}
		
		if(ref.cans[firstNucl] == query.cans[firstNucl] && ref.cans[firstNucl] == -1)
		{
			++TN;
		}
		
	}
	
	
	return std::make_tuple(TP, TN, FP, FN, TPx, TNx, FPx, FNx);
	
}


ConfusionMatrixTuple_blast BLAST::CompareStructures::calcConfusionMatrixNonCans()
{
	using std::cerr;
	using std::endl;
	using std::cout;
	
	int TP, TN, FP, FN, TPx, TNx, FPx, FNx;
	
	TP = TN = FP = FN = TPx = TNx = FPx = FNx = 0;
	
	
	cout << "NonCans:\n";
	cout << "\tref:\n";
	findInteractions.showInformationInRange(findInteractions.m_vRanges[0], ref.nonCans, 0, "n");
	cout << "\tquery\n";
	findInteractions.showInformationInRange(findInteractions.m_vRanges[1], query.nonCans, 1, "n");
	
	
	
	for(size_t firstNucl = 1; firstNucl <=  ref.all.size(); ++firstNucl)
	{
		// the first nucleotide in the pairs are the same
		// so we compare only the second nucleotides
		// [example: vRefSS[i] gives the second nuleotide involved in each pair for the reference structure]
		if(ref.nonCans[firstNucl] == query.nonCans[firstNucl] && ref.nonCans[firstNucl] != -1 )
		{
			if(ref.nonCans[firstNucl] != -2 && ref.nonCans[firstNucl] != -3)
			{
				++TP;
			}
			else
			{
				++TPx;
			}
		}
		
		if(ref.nonCans[firstNucl] != query.nonCans[firstNucl])
		{
			if(ref.nonCans[firstNucl] > 0)
			{
				switch(query.nonCans[firstNucl])
				{
					case -1: ++FN; break;
					case -2: [[fallthrough]];
					case -3: ++FPx;
					default: ++FP; break;
				}
			}
			
			if(ref.nonCans[firstNucl] == -1)
			{
				switch(query.nonCans[firstNucl])
				{
					case -2 : ++FPx; break;
					case -3 : ++FPx; break;
					default : ++FP; break;
				}
			}
			
			if(ref.nonCans[firstNucl] == -2)
			{
				switch(query.nonCans[firstNucl])
				{
					case -1 : ++FNx; break;
					default : ++FPx; break;
				}
			}
			
			if(ref.nonCans[firstNucl] == -3)
			{
				switch(query.nonCans[firstNucl])
				{
					case -1 : ++FNx; break;
					default : ++FPx; break;
				}
			}
		}
		
		if(ref.nonCans[firstNucl] == query.nonCans[firstNucl] && ref.nonCans[firstNucl] == -1)
		{
			++TN;
		}
		
	}
	
	return std::make_tuple(TP, TN, FP, FN, TPx, TNx, FPx, FNx);
	
}


std::vector<ConfusionMatrixTuple_blast> BLAST::CompareStructures::calcConfusionMatrix(fs::path blastPath)
{
	std::vector<ConfusionMatrixTuple_blast> vcmt;
	
	readFile(blastPath);
	
	auto start = requestedInteractions.begin();
	auto end = requestedInteractions.end();
	
	auto itCan = std::find_if(start, end, [](char ch){return (ch == 'c' || ch == 'C');});
	auto itE = std::find_if(start, end, [](char ch){return (ch == 'E' || ch == 'e');});
	if(itCan != end || itE != end)
	{
		vcmt.emplace_back(calcConfusionMatrixCans());
		m_vMax_n_positives.push_back(find_max_n_positives(ref.cans, query.cans));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0, 0, 0, 0, 0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	auto itNonCan = std::find_if(start, end, [](char ch){return (ch == 'n' || ch == 'N');});
	if(itNonCan != end)
	{
		vcmt.emplace_back(calcConfusionMatrixNonCans());
		m_vMax_n_positives.push_back(find_max_n_positives(ref.nonCans, query.nonCans));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0, 0, 0, 0, 0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	auto itAll = std::find_if(start, end, [](char ch){return (ch == 'a'|| ch == 'A');});
	if(itAll != end)
	{
		vcmt.emplace_back(calcConfusionMatrix());
		m_vMax_n_positives.push_back(find_max_n_positives(ref.all, query.all));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0, 0, 0, 0, 0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	return vcmt;
}

// calcSimilarityScores
// using nodiscard attribute to make sure the return value is not being ignored.
std::map<std::string, double> BLAST::CompareStructures::calcSimilarityScores(ConfusionMatrixTuple_blast cmt, int max_n_positives)
{
	std::map<std::string, double> simScoresMap;
	auto [TP, TN, FP, FN, TPx, TNx, FPx, FNx] = cmt;
	// make a object of the SimilarityScores class with the variables of the confiusion matrix tuple
	SimilarityScores simScores(TP + TPx, TN, FP + FPx, FN + FNx);
	
	simScoresMap.insert(std::make_pair("MCC", simScores.calcMCC()));
	simScoresMap.insert(std::make_pair("FScore", simScores.calcFscore()));
	simScoresMap.insert(std::make_pair("FMIndex", simScores.calcFMIndex()));
	simScoresMap.insert(std::make_pair("JIndex", simScores.calcJIndex()));
	simScoresMap.insert(std::make_pair("Precision", simScores.calcPrecision()));
	simScoresMap.insert(std::make_pair("Recall", simScores.calcRecall()));
	simScoresMap.insert(std::make_pair("TP", TP));
	simScoresMap.insert(std::make_pair("TN", TN));
	simScoresMap.insert(std::make_pair("FP", FP));
	simScoresMap.insert(std::make_pair("FN", FN));
	simScoresMap.insert(std::make_pair("TPx", TPx));
	simScoresMap.insert(std::make_pair("TNx", TNx));
	simScoresMap.insert(std::make_pair("FPx", FPx));
	simScoresMap.insert(std::make_pair("FNx", FNx));
	simScoresMap.insert(std::make_pair("Specificity", simScores.calcSpecificty()));
	simScoresMap.insert(std::make_pair("BA", simScores.calcBA()));
	simScoresMap.insert(std::make_pair("FOR", simScores.calcFOR()));
	simScoresMap.insert(std::make_pair("PT", simScores.calcPT()));
	simScoresMap.insert(std::make_pair("CSI", simScores.calcCSI()));
	simScoresMap.insert(std::make_pair("MK", simScores.calcMK()));
	simScoresMap.insert(std::make_pair("JBIndex", static_cast<double>(TP/(max_n_positives + 0.00005))));
	
	
	return simScoresMap;
}


//void BLAST::CompareStructures::writeScores(ConfusionMatrixTuple_blast cmt, fs::path outputPath)
//{
	//std::ofstream outFile(outputPath.string().c_str());
	//auto [TP, TN, FP, FN, TPx, TNx, FPx, FNx] = cmt;
	//std::map<std::string, double> simScoresMap = calcSimilarityScores(cmt);
	//outFile << "Cans,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-" << endl;
	//outFile << "NonCans,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-" << endl;
	//outFile << "All," << simScoresMap["TP"]
			//<< "," << simScoresMap["TN"] << "," << simScoresMap["FP"]
			//<< "," << simScoresMap["FN"]  << "," << TPx
			//<< "," << TNx << "," << FPx << "," << FNx << "," << TP+TN+FP+FN+TPx+TNx+FPx+FNx
			//<< "," << simScoresMap["MCC"] << "," << simScoresMap["FScore"]
			//<< "," << simScoresMap["FMIndex"] << "," << simScoresMap["JIndex"]
			//<< "," << simScoresMap["Precision"] << "," << simScoresMap["Recall"]
			//<< endl;
	
	
	//outFile.close();
	//cout << "your reuslts are in: "  << outputPath.string() << endl;
//}

void BLAST::CompareStructures::writeScores(std::vector<ConfusionMatrixTuple_blast> vcmt, fs::path outputPath)
{
	std::ofstream outFile(outputPath.string().c_str());
	// writing the header to the file
	
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "FSCORE," << "FM_INDEX," << "J_INDEX," 
		    << "PRECISION," << "RECALL," << "Specificity," << "BA," << "FOR,"
		    << "PT," << "CSI," << "MK," << "JBIndex" << std::endl;
	auto itE = std::find_if(requestedInteractions.begin(), requestedInteractions.end(), [](char ch){return (ch == 'E' || ch == 'e');});
	for(size_t i{0}; i < vcmt.size(); ++i)
	{
		auto [TP, TN, FP, FN, TPx, TNx, FPx, FNx] = vcmt[i];
		if(TP == 0 && TN == 0 && FP == 0 && FN == 0 && TPx == 0 && FPx == 0 && FNx == 0)
		{
			continue;
		}
		
		std::map<std::string, double> simScoresMap = calcSimilarityScores(vcmt[i], m_vMax_n_positives[i]);
		
		if(itE == requestedInteractions.end())
		{ // Cans
			
			outFile << vInteractionTypes1[i] << simScoresMap["TP"]
					<< "," << simScoresMap["TN"] << "," << simScoresMap["FP"]
					<< "," << simScoresMap["FN"]  << "," << TP+TN+FP+FN;
					if(simScoresMap["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["MCC"];

					if(simScoresMap["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["FScore"];

					if(simScoresMap["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["FMIndex"];

					if(simScoresMap["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["JIndex"];

					if(simScoresMap["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["Precision"];

					if(simScoresMap["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["Recall"];

					if(simScoresMap["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["Specificity"];

					if(simScoresMap["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["BA"];


					if(simScoresMap["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["FOR"];

					if(simScoresMap["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["PT"];

					if(simScoresMap["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["CSI"];

					if(simScoresMap["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["MK"];

					if(simScoresMap["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["JBIndex"] << endl;
					//<< "," << simScoresMap["MCC"] << "," << simScoresMap["FScore"]
					//<< "," << simScoresMap["FMIndex"] << "," << simScoresMap["JIndex"]
					//<< "," << simScoresMap["Precision"] << "," << simScoresMap["Recall"]
					//<< "," << simScoresMap["Specificity"] << "," << simScoresMap["BA"]
					//<< "," << simScoresMap["FOR"] << "," << simScoresMap["PT"]
					//<< "," << simScoresMap["CSI"] << "," << simScoresMap["MK"]
					//<< "," << simScoresMap["JBIndex"] << endl;
		}
		else
		{ // Cans + Wobble
			//stdExt::tuplePrint(vcmt[i]);
			outFile << vInteractionTypes2[i] << simScoresMap["TP"]
					<< "," << simScoresMap["TN"] << "," << simScoresMap["FP"]
					<< "," << simScoresMap["FN"]  << "," << TP+TN+FP+FN;
					if(simScoresMap["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["MCC"];

					if(simScoresMap["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["FScore"];

					if(simScoresMap["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["FMIndex"];

					if(simScoresMap["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["JIndex"];

					if(simScoresMap["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["Precision"];

					if(simScoresMap["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["Recall"];

					if(simScoresMap["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["Specificity"];

					if(simScoresMap["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["BA"];


					if(simScoresMap["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["FOR"];

					if(simScoresMap["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["PT"];

					if(simScoresMap["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["CSI"];

					if(simScoresMap["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["MK"];

					if(simScoresMap["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMap["JBIndex"] << endl;
					//<< "," << simScoresMap["MCC"] << "," << simScoresMap["FScore"]
					//<< "," << simScoresMap["FMIndex"] << "," << simScoresMap["JIndex"]
					//<< "," << simScoresMap["Precision"] << "," << simScoresMap["Recall"]
					//<< "," << simScoresMap["Specificity"] << "," << simScoresMap["BA"]
					//<< "," << simScoresMap["FOR"] << "," << simScoresMap["PT"]
					//<< "," << simScoresMap["CSI"] << "," << simScoresMap["MK"]
					//<< "," << simScoresMap["JBIndex"] << endl;
		}
	}
	
	
	outFile.close();
	cout << "your reuslts are in: "  << outputPath.string() << endl;
}


	
// a template helper function and functions for printing tuples
template<typename Tuple, size_t... Indices>
void tuplePrint(const Tuple& t, std::index_sequence<Indices...>)
{
	((cout << std::get<Indices>(t) << " "), ...);
}

template<typename... Args>
void tuplePrint(const std::tuple<Args...>& t)
{
	tuplePrint(t, std::index_sequence_for<Args...>());
	cout << endl;
}

//#define _MAIN_
#ifdef _MAIN_
using namespace std;
int main()
{
	BLAST::CompareStructures cs;
	
	//cs.findInteractions.set_isWobble_canonical(true);
	//cs.readFile("blast_example.txt");
	cs.requestedInteractions = "cn";
	auto vcmt = cs.calcConfusionMatrix("blast_example.txt");
	cs.writeScores(vcmt, "FINAL_results.csv");
	//auto confusionMatrix = cs.calcConfusionMatrix(/*"blast_example.txt"*/);
	////tuplePrint(confusionMatrix);
	
	//auto confusionMatrixCans = cs.calcConfusionMatrixCans(/*"blast_example.txt"*/);
	////tuplePrint(confusionMatrixCans);
	
	//auto confusionMatrixNonCans = cs.calcConfusionMatrixNonCans(/*"blast_example.txt"*/);
	//tuplePrint(confusionMatrixNonCans);
	cout << endl;
}

#endif //_MAIN_
