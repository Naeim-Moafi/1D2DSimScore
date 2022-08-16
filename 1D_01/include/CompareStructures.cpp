#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

using std::endl;
const std::vector<std::string> vInteractionTypes1 = {"Cans,", "wobbles,", "NonCans,", "AllBasePairs,"};
const std::vector<std::string> vInteractionTypes2 = {"Cans + wobbles,", "wobbles,", "NonCans,", "AllBasePairs,"};

//constructors
void Binary::CompareStructures::readInputFiles(const fs::path& refPath, const fs::path& queryPath)
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
	
	
	// checking if the number of the chains for the inputs are the same
	if(findInteractionRef.m_vChains_binary.size() != findInteractionQuery.m_vChains_binary.size())
	{
		throw std::invalid_argument("The size of the secondary structures are not the same");
	}
	
	// checking if the size of the corresponding chain in both inputs are the same
	for(size_t i_chain { 0 }; i_chain < findInteractionRef.m_vChains_binary.size(); ++i_chain)
	{
		if(findInteractionRef.m_vChains_binary[i_chain].size() != findInteractionQuery.m_vChains_binary[i_chain].size())
		{
			std::ostringstream oss_err_message;
			oss_err_message <<"the size of the chain " << i_chain << " in reference and query is different\nPlease check your inputs\n";
			throw std::invalid_argument(oss_err_message.str());
		}
	}
}

ConfusionMatrixTuple Binary::CompareStructures::calcConfusionMatrix()
{	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;

	for(size_t i_res { 0 }; i_res < findInteractionRef.binaryInterface.size(); ++i_res)
	{
		std::string ref(1,findInteractionRef.binaryInterface[i_res]);
		std::string query(1,findInteractionQuery.binaryInterface[i_res]);
		if(ref == query && (ref == "X" || ref == "1"))
		{
			++findInteractionRef.m_number_positives;
			++findInteractionQuery.m_number_positives;
			++TP;
		}
		
		if(ref != query)
		{
			if(ref == "X" || ref == "1")
			{
				++findInteractionRef.m_number_positives;
				++FN;
			}
			else
			{
				++findInteractionQuery.m_number_positives;
				++FP;
			}
		}
		
		
		if(ref == query && (ref == "." || ref == "0"))
		{
			++TN;
		}
	}
	
	if(findInteractionRef.m_number_positives > findInteractionQuery.m_number_positives)
	{
		m_max_n_positives = findInteractionRef.m_number_positives;
	}
	else
	{
		m_max_n_positives = findInteractionQuery.m_number_positives;
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

std::map<std::string, double> Binary::CompareStructures::calcSimilarityScores(ConfusionMatrixTuple cmt)
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
	simScoresMap.insert(std::make_pair("JBIndex", static_cast<double>(TP/(m_max_n_positives + 0.00005))));

	
	return simScoresMap;
}


void Binary::CompareStructures::writeScores(ConfusionMatrixTuple cmt, fs::path outputPath)
{
	std::ofstream outFile(outputPath.string().c_str());
	auto [TP, TN, FP, FN] = cmt;
	std::map<std::string, double> simScoresMap = calcSimilarityScores(cmt);
	
	
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "F1," << "FM," << "J," 
		    << " PRECISION," << " RECALL," << " Specificity," << "BA," << "FOR,"
		    << " PT," << "CSI," << "MK," << "B" << std::endl;
	outFile << "Binary," << simScoresMap["TP"]
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
	
	
	outFile.close();
	cout << "your reuslts are in: "  << outputPath.string() << endl;
}

