#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

using std::endl;

// readFiles
void CMO::CompareStructures::readFiles(const fs::path& refPath, const fs::path& queryPath)
{
	try
	{
		findInteractionRef.readInputFile(refPath);
		findInteractionQuery.readInputFile(queryPath);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
		exit(EXIT_FAILURE);
	}
	
	
	//checking if the size of the structures are the same
	if(findInteractionRef.m_Matrix.size() != findInteractionQuery.m_Matrix.size())
	{
		throw std::invalid_argument("The size of the secondary structures are not the same");
	}
}

void CMO::CompareStructures::find_max_n_positives(const Matrix& ref_, const Matrix& query_)
{
	int np_ref = 0;
	int np_query = 0;
	for(size_t i { 0 }; i < ref_.size(); ++i)
	{
		for(size_t j { 0 }; j < ref_.size(); ++j)
		{
			if(ref_[i][j] != '0')
			{
				++np_ref;
			}
		}
	}
	
	for(size_t i { 0 }; i < ref_.size(); ++i)
	{
		for(size_t j { 0 }; j < ref_.size(); ++j)
		{
			if(query_[i][j] != '0')
			{
				++np_query;
			}
		}
	}
	
	
	if(np_ref > np_query)
	{
		m_max_n_positives = np_ref;
	}
	
	m_max_n_positives = np_query;
}


ConfusionMatrixTuple CMO::CompareStructures::calcConfusionMatrix()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	int length = findInteractionRef.m_Matrix.size();
	
	
	for(size_t firstNucl = 0;  firstNucl <  findInteractionRef.m_Matrix.size(); ++firstNucl)
	{
		for(size_t secondNucle = firstNucl + 1; secondNucle < findInteractionRef.m_Matrix.size(); ++secondNucle)
		{
			if(findInteractionRef.m_Matrix[firstNucl][secondNucle] == findInteractionQuery.m_Matrix[firstNucl][secondNucle] && findInteractionRef.m_Matrix[firstNucl][secondNucle] != '0')
			{
				TP++;
			}
			if(findInteractionRef.m_Matrix[firstNucl][secondNucle] != findInteractionQuery.m_Matrix[firstNucl][secondNucle])
			{
				if(findInteractionRef.m_Matrix[firstNucl][secondNucle] == '0')
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
			//if(findInteractionRef.m_Matrix[firstNucl][secondNucle] == findInteractionQuery.m_Matrix[firstNucl][secondNucle] && findInteractionRef.m_Matrix[firstNucl][secondNucle] == '0')
			//{
				//++TN;
			//}
		}
	}
	int possible_intercations_number = (length - 1) * (length - 2) / 2;
	TN = possible_intercations_number - (TP + FP + FN);
	return std::make_tuple(TP, TN, FP, FN);	
}

std::map<std::string, double> CMO::CompareStructures::calcSimilarityScores(ConfusionMatrixTuple cmt, int max_n_positives)
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
	simScoresMap.insert(std::make_pair("Specificity", simScores.calcSpecificty()));
	simScoresMap.insert(std::make_pair("BA", simScores.calcBA()));
	simScoresMap.insert(std::make_pair("FOR", simScores.calcFOR()));
	simScoresMap.insert(std::make_pair("PT", simScores.calcPT()));
	simScoresMap.insert(std::make_pair("CSI", simScores.calcCSI()));
	simScoresMap.insert(std::make_pair("MK", simScores.calcMK()));
	//simScoresMap.insert(std::make_pair("JBIndex", static_cast<double>(TP/(max_n_positives + 0.00005))));

	
	return simScoresMap;
}




//void SS::CompareStructures::writeScores(ConfusionMatrixTuple cmt, fs::path outputPath)
//{
	//std::ofstream outFile(outputPath.string().c_str());
	//auto [TP, TN, FP, FN] = cmt;
	//std::map<std::string, double> simScoresMap = calcSimilarityScores(cmt);
	
	//outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	//<< " MCC,"  << "FSCORE," << "FM_INDEX," << "J_INDEX," 
		    //<< " PRECISION," << " RECALL"<< std::endl;
	//outFile << "Cans,-,-,-,-,-,-,-,-,-,-,-" << endl;
	//outFile << "NonCans,-,-,-,-,-,-,-,-,-,-,-" << endl;
	//outFile << "AllBasePairs," << simScoresMap["TP"]
			//<< "," << simScoresMap["TN"] << "," << simScoresMap["FP"]
			//<< "," << simScoresMap["FN"]  << "," << TP+TN+FP+FN
			//<< "," << simScoresMap["MCC"] << "," << simScoresMap["FScore"]
			//<< "," << simScoresMap["FMIndex"] << "," << simScoresMap["JIndex"]
			//<< "," << simScoresMap["Precision"] << "," << simScoresMap["Recall"]
			//<< endl;
	
	//outFile.close();
	//cout << "your reuslts are in: "  << outputPath.string() << endl;
//}



void CMO::CompareStructures::writeScores(ConfusionMatrixTuple cmt, fs::path outputPath)
{
	std::ofstream outFile(outputPath.string().c_str());
	// writing the header to the file
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "F1," << "FM," << "J," 
		    << " PRECISION," << " RECALL," << " Specificity," << " BA," << " FOR,"
		    << "PT," << "CSI," << "MK,"/* << "B"*/ << std::endl;

	auto [TP, TN, FP, FN] = cmt;
	std::map<std::string, double> simScoresMap = calcSimilarityScores(cmt, m_max_n_positives);
	outFile << "All," << simScoresMap["TP"]
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

			//if(simScoresMap["JBIndex"] == -1.1) outFile << ",-";
			//else outFile << "," << simScoresMap["JBIndex"] << endl;
	
	outFile.close();
	cout << "your reuslts are in: "  << outputPath.string() << endl;
}


//#define _MAIN_
#ifdef _MAIN_
using namespace std;
int main()
{
	CMO::CompareStructures cs;
	cs.readFiles("1e6v_1_All.map","1e6v_1_All.map");
	auto cmt = cs.calcConfusionMatrix();
	cs.writeScores(cmt, "Final_resultes.csv");
}

#endif //_MAIN_
