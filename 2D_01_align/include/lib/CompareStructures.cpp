#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

using std::endl;
const std::vector<std::string> vInteractionTypes = {"Can,", "NonCan,", "AllBasePairs,"};

//constructors
SS::CompareStructures::CompareStructures()
{
	findInteractionRef.set_withSeq(false);
	findInteractionRef.set_isWobble_canonical(false);
	findInteractionQuery.set_withSeq(false);
	findInteractionQuery.set_isWobble_canonical(false);
	requestedInteractions = "A";
}

SS::CompareStructures::CompareStructures(bool withSeq, bool isWobble_canonical) 
{
	findInteractionRef.set_withSeq(withSeq);
	findInteractionRef.set_isWobble_canonical(isWobble_canonical);
	findInteractionQuery.set_withSeq(withSeq);
	findInteractionQuery.set_isWobble_canonical(isWobble_canonical);
}

// readStructuresSS
void SS::CompareStructures::readStructures(fs::path refPath, fs::path queryPath)
{
	try
	{
		vRef = findInteractionRef.readInputFile(refPath);
		vQuery = findInteractionQuery.readInputFile(queryPath);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	
	//checking if the size of the structures are the same
	if(findInteractionRef.ssSize != findInteractionQuery.ssSize)
	{
		throw std::invalid_argument("The size of the secondary structures are not the same");
	}
}

// calcConfusionMatrixFromSS
// using nodiscard attribute to make sure the return value is not being ignored. 
ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrix(fs::path refPath, fs::path queryPath)
{
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	//reading ss files
	try
	{
		readStructures(refPath, queryPath);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	
	for(size_t firstNucl = 1; firstNucl <=  vRef.size(); ++firstNucl)
	{
		// the first nucleotide in the pairs are the same
		// so we compare only the second nucleotides
		// [example: vRefSS[i] gives the second nuleotide involved in each pair for the reference structure]
		if(vRef[firstNucl] == vQuery[firstNucl] && vRef[firstNucl] != -1)
		{
			++TP;
		}
		
		//if(vRefSS[firstNucl] != vQuerySS[firstNucl] && vRefSS[firstNucl] == -1)
		//{
			//++FP;
		//}
		
		//if(vRefSS[firstNucl] != vQuerySS[firstNucl] && vRefSS[firstNucl] != -1)
		//{
			//++FN;
		//}
		
		if(vRef[firstNucl] != vQuery[firstNucl] /*&& vRefSS[firstNucl]*/)
		{
			if(vRef[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(vQuery[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
		
		if(vRef[firstNucl] == vQuery[firstNucl] && vRef[firstNucl] == -1)
		{
			++TN;
		}
		
	}
	
	return std::make_tuple(TP, TN, FP, FN);
	
}

std::vector<ConfusionMatrixTuple> SS::CompareStructures::calcConfusionMatrix(fs::path refPath, fs::path queryPath, fs::path seqPath)
{
	std::vector<ConfusionMatrixTuple> vcmt;
	
	//reading ss files
	try
	{
		readStructures(refPath, queryPath);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	sequence = findInteractionRef.readSeqFile(seqPath);
	findInteractionRef.mapSS2seq(sequence, vRef);
	findInteractionQuery.mapSS2seq(sequence, vQuery);
	
	auto start = requestedInteractions.begin();
	auto end = requestedInteractions.end();
	
	auto itCan = std::find_if(start, end, [](char ch){return (ch == 'c' || ch == 'C');});
	if(itCan != end)
	{
		vcmt.emplace_back(calcConfusionMatrixCans());
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0, 0, 0, 0));
	}
	
	auto itNonCan = std::find_if(start, end, [](char ch){return (ch == 'n' || ch == 'N');});
	if(itNonCan != end)
	{
		vcmt.emplace_back(calcConfusionMatrixNonCans());
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0, 0, 0, 0));
	}
	
	auto itAll = std::find_if(start, end, [](char ch){return (ch == 'a'|| ch == 'A');});
	if(itAll != end)
	{
		vcmt.emplace_back(calcConfusionMatrix(refPath, queryPath));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0, 0, 0, 0));
	}
	
	return vcmt;
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixCans()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	
	for(size_t firstNucl = 1; firstNucl <=  findInteractionRef.ssMap_can.size(); ++firstNucl)
	{
		if(findInteractionRef.ssMap_can[firstNucl] == findInteractionQuery.ssMap_can[firstNucl] && findInteractionRef.ssMap_can[firstNucl] != -1)
		{
			TP++;
		}
		if(findInteractionRef.ssMap_can[firstNucl] != findInteractionQuery.ssMap_can[firstNucl] /*&& vRefSS[firstNucl]*/)
		{
			if(findInteractionRef.ssMap_can[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(findInteractionQuery.ssMap_can[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
		
		if(findInteractionRef.ssMap_can[firstNucl] == findInteractionQuery.ssMap_can[firstNucl] && findInteractionRef.ssMap_can[firstNucl] == -1)
		{
			++TN;
		}
	}
		
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixNonCans()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	
	
	for(size_t firstNucl = 1; firstNucl <=  findInteractionRef.ssMap_noncan.size(); ++firstNucl)
	{
		if(findInteractionRef.ssMap_noncan[firstNucl] == findInteractionQuery.ssMap_noncan[firstNucl] && findInteractionRef.ssMap_noncan[firstNucl] != -1)
		{
			TP++;
		}
		if(findInteractionRef.ssMap_noncan[firstNucl] != findInteractionQuery.ssMap_noncan[firstNucl] /*&& vRefSS[firstNucl]*/)
		{
			if(findInteractionRef.ssMap_noncan[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(findInteractionQuery.ssMap_noncan[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
		
		if(findInteractionRef.ssMap_noncan[firstNucl] == findInteractionQuery.ssMap_noncan[firstNucl] && findInteractionRef.ssMap_noncan[firstNucl] == -1)
		{
			++TN;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

// calcSimilarityScores
// using nodiscard attribute to make sure the return value is not being ignored.
[[nodiscard]]std::map<std::string, double> SS::CompareStructures::calcSimilarityScores(ConfusionMatrixTuple cmt)
{
	std::map<std::string, double> simScoresMap;
	// make a object of the SimilarityScores class with the variables of the confiusion matrix tuple
	auto [TP, TN, FP, FN] = cmt;
	
	SimilarityScores simScores(TP, TN, FP, FN);
	
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
	
	return simScoresMap;
}


void SS::CompareStructures::writeScores(ConfusionMatrixTuple cmt, fs::path outputPath)
{
	std::ofstream outFile(outputPath.string().c_str());
	auto [TP, TN, FP, FN] = cmt;
	std::map<std::string, double> simScoresMap = calcSimilarityScores(cmt);
	
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "FSCORE," << "FM_INDEX," << "J_INDEX," 
		    << " PRECISION," << " RECALL"<< std::endl;
	outFile << "Cans,-,-,-,-,-,-,-,-,-,-,-" << endl;
	outFile << "NonCans,-,-,-,-,-,-,-,-,-,-,-" << endl;
	outFile << "AllBasePairs," << simScoresMap["TP"]
			<< "," << simScoresMap["TN"] << "," << simScoresMap["FP"]
			<< "," << simScoresMap["FN"]  << "," << TP+TN+FP+FN
			<< "," << simScoresMap["MCC"] << "," << simScoresMap["FScore"]
			<< "," << simScoresMap["FMIndex"] << "," << simScoresMap["JIndex"]
			<< "," << simScoresMap["Precision"] << "," << simScoresMap["Recall"]
			<< endl;
	
	
	outFile.close();
	cout << "your reuslts are in: "  << outputPath.string() << endl;
}

void SS::CompareStructures::writeScores(std::vector<ConfusionMatrixTuple> vcmt, fs::path outputPath)
{
	std::ofstream outFile(outputPath.string().c_str());
	// writing the header to the file
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "FSCORE," << "FM_INDEX," << "J_INDEX," 
		    << " PRECISION," << " RECALL"<< std::endl;
	
	for(size_t i{0}; i < vcmt.size(); ++i)
	{
		auto [TP, TN, FP, FN] = vcmt[i];
		if(TP == 0 && TN == 0 && FP == 0 && FN == 0)
		{
			outFile << vInteractionTypes[i] << "-,-,-,-,-,-,-,-,-,-,-" << endl;
		}
		else
		{
			std::map<std::string, double> simScoresMap = calcSimilarityScores(vcmt[i]);
			outFile << vInteractionTypes[i] << simScoresMap["TP"]
					<< "," << simScoresMap["TN"] << "," << simScoresMap["FP"]
					<< "," << simScoresMap["FN"]  << "," << TP+TN+FP+FN
					<< "," << simScoresMap["MCC"] << "," << simScoresMap["FScore"]
					<< "," << simScoresMap["FMIndex"] << "," << simScoresMap["JIndex"]
					<< "," << simScoresMap["Precision"] << "," << simScoresMap["Recall"]
					<< endl;
		}
	}
	outFile.close();
	cout << "your reuslts are in: "  << outputPath.string() << endl;
}


//#define _MAIN_
#ifdef _MAIN_
using namespace std;
int main()
{
	SS::CompareStructures cs;
	
	//cs.findInteractionRef.set_isWobble_canonical(false);
	//cs.findInteractionQuery.set_isWobble_canonical(false);
	
	//write the SSmaps of the ref and query
	vector<ConfusionMatrixTuple> vcmt = cs.calcConfusionMatrix("dotBracketRef.SS", "dotBracketQuery.SS", "SeqForDotBracket.seq");
	cs.writeScores(vcmt, "Final_resultes.csv");
	
	//cout << "TP: " << get<0>(cmtSS) << " TN: " << get<1>(cmtSS) << " FP: " << get<2>(cmtSS) << " FN: " << get<3>(cmtSS) << endl;
	
	//SimilarityScores sscSS(get<0>(cmtSS), get<1>(cmtSS), get<2>(cmtSS), get<3>(cmtSS));
	//cs.writeScores(vcmtClaRNA, "finalOut.csv");
	//ofstream outFile("finalOut.csv");
	//cout << "ALL:\n";
	
	//map<string, double> simScoresMapA = cs.calcSimilarityScores(vcmtClaRNA[0]);
	//cout << "TP: " << " TN: " << " FP: " << simScoresMapA["FP"] << " FN: " 
		 //<< " MCC: "  << " FSCORE: " << " FM INDESX: " << " J INDEX: " 
		 //<< " PRECISION: " << " RECALL: "<< endl;
		 
	//cout << simScoresMapA["TP"] << " " << simScoresMapA["TN"] << " " << simScoresMapA["FP"] << " " 
		 //<< simScoresMapA["FN"] << " " << simScoresMapA["MCC"] << " " << simScoresMapA["FScore"] << " " 
		 //<< simScoresMapA["FMIndex"] << " " << simScoresMapA["JIndex"] << " " << simScoresMapA["Precision"] << " " 
		 //<< simScoresMapA["Recall"] << endl;
		 

	
	//outFile << "All," << simScoresMapA["TP"] << "," << simScoresMapA["TN"] << "," << simScoresMapA["FP"] << "," 
		 //<< simScoresMapA["FN"] << "," << simScoresMapA["MCC"] << "," << simScoresMapA["FScore"] << "," 
		 //<< simScoresMapA["FMIndex"] << "," << simScoresMapA["JIndex"] << "," << simScoresMapA["Precision"] << "," 
		 //<< simScoresMapA["Recall"] << endl;
	
	
	
	
	//cout << "\nCanonicals:\n";
	
	//map<string, double> simScoresMapC = cs.calcSimilarityScores(vcmtClaRNA[1]);
	//cout << simScoresMapC["TP"] << " " << simScoresMapC["TN"] << " " << simScoresMapC["FP"] << " " 
		 //<< simScoresMapC["FN"] << " " << simScoresMapC["MCC"] << " " << simScoresMapC["FScore"] << " " 
		 //<< simScoresMapC["FMIndex"] << " " << simScoresMapC["JIndex"] << " " << simScoresMapC["Precision"] << " " 
		 //<< simScoresMapC["Recall"] << endl;
	
	
	//outFile << "Cans," << simScoresMapC["TP"] << "," << simScoresMapC["TN"] << "," << simScoresMapC["FP"] << "," 
		 //<< simScoresMapC["FN"] << "," << simScoresMapC["MCC"] << "," << simScoresMapC["FScore"] << "," 
		 //<< simScoresMapC["FMIndex"] << "," << simScoresMapC["JIndex"] << "," << simScoresMapC["Precision"] << "," 
		 //<< simScoresMapC["Recall"] << endl;
	
	//cout << "\nNoncanonicals:\n";
	
	//map<string, double> simScoresMapN = cs.calcSimilarityScores(vcmtClaRNA[2]);
	//cout << simScoresMapN["TP"] << " " << simScoresMapN["TN"] << " " << simScoresMapN["FP"] << " " 
		 //<< simScoresMapN["FN"] << " " << simScoresMapN["MCC"] << " " << simScoresMapN["FScore"] << " "
		 //<< simScoresMapN["FMIndex"] << " " << simScoresMapN["JIndex"] << " " << simScoresMapN["Precision"] << " " 
		 //<< simScoresMapN["Recall"] << endl;
	
	
	//outFile << "NonCans," << simScoresMapN["TP"] << "," << simScoresMapN["TN"] << "," << simScoresMapN["FP"] << "," 
		 //<< simScoresMapN["FN"] << "," << simScoresMapN["MCC"] << "," << simScoresMapN["FScore"] << "," 
		 //<< simScoresMapN["FMIndex"] << "," << simScoresMapN["JIndex"] << "," << simScoresMapN["Precision"] << "," 
		 //<< simScoresMapN["Recall"] << endl;
	
	//cout << "\nStacks:\n";
	
	//map<string, double> simScoresMapS = cs.calcSimilarityScores(vcmtClaRNA[3]);
	//cout << simScoresMapS["TP"] << " " << simScoresMapS["TN"] << " " << simScoresMapS["FP"] << " " 
		 //<< simScoresMapS["FN"] << " " << simScoresMapS["MCC"] << " " << simScoresMapS["FScore"] << " " 
		 //<< simScoresMapS["FMIndex"] << " " << simScoresMapS["JIndex"] << " " << simScoresMapS["Precision"] << " " 
		 //<< simScoresMapS["Recall"] << endl;
	
	//outFile << "stacks," << simScoresMapS["TP"] << "," << simScoresMapS["TN"] << "," << simScoresMapS["FP"] << "," 
		 //<< simScoresMapS["FN"] << "," << simScoresMapS["MCC"] << "," << simScoresMapS["FScore"] << "," 
		 //<< simScoresMapS["FMIndex"] << "," << simScoresMapS["JIndex"] << "," << simScoresMapS["Precision"] << "," 
		 //<< simScoresMapS["Recall"] << endl;
	
	//outFile.close();
}

#endif //_MAIN_
