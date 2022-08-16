#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

using std::endl;
const std::vector<std::string> vInteractionTypes1 = {"Cans,", "wobbles,", "NonCans,", "AllBasePairs,"};
const std::vector<std::string> vInteractionTypes2 = {"Cans + wobbles,", "wobbles,", "NonCans,", "AllBasePairs,"};

//constructors
SS::CompareStructures::CompareStructures()
{
	findInteractionRef.set_isWobble_canonical(false);
	findInteractionQuery.set_isWobble_canonical(false);
	m_requestedInteractions = "A";
}

SS::CompareStructures::CompareStructures(bool isWobble_canonical, bool withSeq) 
{
	findInteractionRef.set_isWobble_canonical(isWobble_canonical);
	findInteractionRef.set_withSeq(withSeq);
	
	findInteractionQuery.set_isWobble_canonical(isWobble_canonical);
	findInteractionQuery.set_withSeq(withSeq);
}

// readStructuresSS
void SS::CompareStructures::readStructures(const fs::path& refPath, const fs::path& queryPath)
{
	try
	{
		findInteractionRef.init_structure(refPath);
		findInteractionQuery.init_structure(queryPath);
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

int SS::CompareStructures::find_max_n_positives(const SSMap& ref_, const SSMap& query_) const
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


int SS::CompareStructures::find_max_n_positives(const SSMatrix& ref_, const SSMatrix& query_) const
{
	int np_ref = 0;
	int np_query = 0;
	for(size_t i { 0 }; i < ref_.size(); ++i)
	{
		for(size_t j { 0 }; j < ref_.size(); ++j)
		{
			if(ref_[i][j] == 1)
			{
				++np_ref;
			}
		}
	}
	
	for(size_t i { 0 }; i < ref_.size(); ++i)
	{
		for(size_t j { 0 }; j < ref_.size(); ++j)
		{
			if(query_[i][j] == 1)
			{
				++np_query;
			}
		}
	}
	
	
	if(np_ref > np_query)
	{
		return np_ref;
	}
	
	return np_query;
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixAllVector()
{	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	for(size_t firstNucl = 1; firstNucl <=  findInteractionRef.m_SSMap.size(); ++firstNucl)
	{
		if(findInteractionRef.m_SSMap[firstNucl] == findInteractionQuery.m_SSMap[firstNucl] && findInteractionRef.m_SSMap[firstNucl] != -1)
		{
			++TP;
		}
		
		if(findInteractionRef.m_SSMap[firstNucl] != findInteractionQuery.m_SSMap[firstNucl])
		{
			if(findInteractionRef.m_SSMap[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(findInteractionQuery.m_SSMap[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
			
		}
		
		if(findInteractionRef.m_SSMap[firstNucl] == findInteractionQuery.m_SSMap[firstNucl] && findInteractionRef.m_SSMap[firstNucl] == -1)
		{
			++TN;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixAllMatrix()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	int length = findInteractionRef.m_Matrix_all.size();
	
	
	for(size_t firstNucl = 0;  firstNucl <  findInteractionRef.m_Matrix_all.size(); ++firstNucl)
	{
		for(size_t secondNucle = firstNucl + 1; secondNucle < findInteractionRef.m_Matrix_all.size(); ++secondNucle)
		{
			if(findInteractionRef.m_Matrix_all[firstNucl][secondNucle] == findInteractionQuery.m_Matrix_all[firstNucl][secondNucle] && findInteractionRef.m_Matrix_all[firstNucl][secondNucle] != 0)
			{
				TP++;
			}
			if(findInteractionRef.m_Matrix_all[firstNucl][secondNucle] != findInteractionQuery.m_Matrix_all[firstNucl][secondNucle])
			{
				if(findInteractionRef.m_Matrix_all[firstNucl][secondNucle] == 0)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
	}
	
	int possible_intercations_number = (length - 1) * (length - 2) / 2;
	TN = possible_intercations_number - (TP + FP + FN);
	return std::make_tuple(TP, TN, FP, FN);	
}


ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixCansVector()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	for(size_t firstNucl = 1; firstNucl <=  findInteractionRef.m_SSMap_can.size(); ++firstNucl)
	{
		if(findInteractionRef.m_SSMap_can[firstNucl] == findInteractionQuery.m_SSMap_can[firstNucl] && findInteractionRef.m_SSMap_can[firstNucl] != -1)
		{
			TP++;
		}
		if(findInteractionRef.m_SSMap_can[firstNucl] != findInteractionQuery.m_SSMap_can[firstNucl] /*&& vRefSS[firstNucl]*/)
		{
			if(findInteractionRef.m_SSMap_can[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(findInteractionQuery.m_SSMap_can[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
		
		if(findInteractionRef.m_SSMap_can[firstNucl] == findInteractionQuery.m_SSMap_can[firstNucl] && findInteractionRef.m_SSMap_can[firstNucl] == -1)
		{
			++TN;
		}
	}
		
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixWobblesVector()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	
	for(size_t firstNucl = 1; firstNucl <=  findInteractionRef.m_SSMap_can.size(); ++firstNucl)
	{
		if(findInteractionRef.m_SSMap_wobble[firstNucl] == findInteractionQuery.m_SSMap_wobble[firstNucl] && findInteractionRef.m_SSMap_wobble[firstNucl] != -1)
		{
			TP++;
		}
		if(findInteractionRef.m_SSMap_wobble[firstNucl] != findInteractionQuery.m_SSMap_wobble[firstNucl] /*&& vRefSS[firstNucl]*/)
		{
			if(findInteractionRef.m_SSMap_wobble[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(findInteractionQuery.m_SSMap_wobble[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
		
		if(findInteractionRef.m_SSMap_wobble[firstNucl] == findInteractionQuery.m_SSMap_wobble[firstNucl] && findInteractionRef.m_SSMap_wobble[firstNucl] == -1)
		{
			++TN;
		}
	}
		
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixCansMatrix()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	int length = findInteractionRef.m_Matrix_can.size();
	
	for(size_t firstNucl = 0;  firstNucl <  findInteractionRef.m_SSMap_can.size(); ++firstNucl)
	{
		for(size_t secondNucle = firstNucl + 1; secondNucle < findInteractionRef.m_Matrix_can.size(); ++secondNucle)
		{
			if(findInteractionRef.m_Matrix_can[firstNucl][secondNucle] == findInteractionQuery.m_Matrix_can[firstNucl][secondNucle] && findInteractionRef.m_Matrix_can[firstNucl][secondNucle] != 0)
			{
				TP++;
			}
			if(findInteractionRef.m_Matrix_can[firstNucl][secondNucle] != findInteractionQuery.m_Matrix_can[firstNucl][secondNucle])
			{
				if(findInteractionRef.m_Matrix_can[firstNucl][secondNucle] == 0)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
	}
	
	int possible_intercations_number = (length - 1) * (length - 2) / 2;
	TN = possible_intercations_number - (TP + FP + FN);
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixWobblesMatrix()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	int length = findInteractionRef.m_Matrix_can.size();
	
	for(size_t firstNucl = 0;  firstNucl <  findInteractionRef.m_SSMap_wobble.size(); ++firstNucl)
	{
		for(size_t secondNucle = firstNucl + 1; secondNucle < findInteractionRef.m_Matrix_can.size(); ++secondNucle)
		{
			if(findInteractionRef.m_Matrix_wobble[firstNucl][secondNucle] == findInteractionQuery.m_Matrix_wobble[firstNucl][secondNucle] && findInteractionRef.m_Matrix_wobble[firstNucl][secondNucle] != 0)
			{
				TP++;
			}
			if(findInteractionRef.m_Matrix_wobble[firstNucl][secondNucle] != findInteractionQuery.m_Matrix_wobble[firstNucl][secondNucle])
			{
				if(findInteractionRef.m_Matrix_wobble[firstNucl][secondNucle] == 0)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
	}
	
	int possible_intercations_number = (length - 1) * (length - 2) / 2;
	TN = possible_intercations_number - (TP + FP + FN);
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixNonCansVector()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	
	
	for(size_t firstNucl = 1; firstNucl <=  findInteractionRef.m_SSMap_noncan.size(); ++firstNucl)
	{
		if(findInteractionRef.m_SSMap_noncan[firstNucl] == findInteractionQuery.m_SSMap_noncan[firstNucl] && findInteractionRef.m_SSMap_noncan[firstNucl] != -1)
		{
			TP++;
		}
		if(findInteractionRef.m_SSMap_noncan[firstNucl] != findInteractionQuery.m_SSMap_noncan[firstNucl] /*&& vRefSS[firstNucl]*/)
		{
			if(findInteractionRef.m_SSMap_noncan[firstNucl] == -1)
			{
				++FP;
			}
			else
			{
				if(findInteractionQuery.m_SSMap_noncan[firstNucl] != -1)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
		
		if(findInteractionRef.m_SSMap_noncan[firstNucl] == findInteractionQuery.m_SSMap_noncan[firstNucl] && findInteractionRef.m_SSMap_noncan[firstNucl] == -1)
		{
			++TN;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

ConfusionMatrixTuple SS::CompareStructures::calcConfusionMatrixNonCansMatrix()
{
	int TP, TN, FP, FN;
	TP = TN = FP = FN = 0;
	int length = findInteractionRef.m_Matrix_noncan.size();
	
	for(size_t firstNucl = 0;  firstNucl <  findInteractionRef.m_SSMap_noncan.size(); ++firstNucl)
	{
		for(size_t secondNucle = firstNucl + 1; secondNucle < findInteractionRef.m_Matrix_noncan.size(); ++secondNucle)
		{
			if(findInteractionRef.m_Matrix_noncan[firstNucl][secondNucle] == findInteractionQuery.m_Matrix_noncan[firstNucl][secondNucle] && findInteractionRef.m_Matrix_noncan[firstNucl][secondNucle] != 0)
			{
				TP++;
			}
			if(findInteractionRef.m_Matrix_noncan[firstNucl][secondNucle] != findInteractionQuery.m_Matrix_noncan[firstNucl][secondNucle])
			{
				if(findInteractionRef.m_Matrix_noncan[firstNucl][secondNucle] == 0)
				{
					++FP;
				}
				else
				{
					++FN;
				}
			}
		}
	}
	
	int possible_intercations_number = (length - 1) * (length - 2) / 2;
	TN = possible_intercations_number - (TP + FP + FN);
	return std::make_tuple(TP, TN, FP, FN);
}

std::vector<ConfusionMatrixTuple> SS::CompareStructures::calcConfusionMatrixVector()
{
	auto start = m_requestedInteractions.begin();
	auto end=  m_requestedInteractions.end();
	std::vector<ConfusionMatrixTuple> vcmt;
	
	// check the requested scores
	//Cans
	auto itE = std::find_if(start, end, [](char ch){return (ch == 'E' || ch == 'e');});
	auto itC = std::find_if(start, end, [](char ch){return (ch == 'C' || ch == 'c');});
	if(itC != end || itE != end)
	{
		vcmt.emplace_back(calcConfusionMatrixCansVector());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_SSMap_can, findInteractionQuery.m_SSMap_can));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
		
	//Wobbles
	auto itW = std::find_if(start, end, [](char ch){return (ch == 'W' || ch == 'w');});
	if(itW != end)
	{
		vcmt.emplace_back(calcConfusionMatrixWobblesVector());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_SSMap_wobble, findInteractionQuery.m_SSMap_wobble));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
	
	//NonCans
	auto itN = std::find_if(start, end, [](char ch){return (ch == 'N' || ch == 'n');});
	if(itN != end)
	{
		vcmt.emplace_back(calcConfusionMatrixNonCansVector());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_SSMap_noncan, findInteractionQuery.m_SSMap_noncan));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
	
	
	//All
	auto itA = std::find_if(start, end, [](char ch){return (ch == 'A' || ch == 'a');});
	if(itA != end)
	{
		vcmt.emplace_back(calcConfusionMatrixAllVector());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_SSMap, findInteractionQuery.m_SSMap));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
	
	return vcmt;
}


std::vector<ConfusionMatrixTuple> SS::CompareStructures::calcConfusionMatrixMatrix()
{
	auto start = m_requestedInteractions.begin();
	auto end=  m_requestedInteractions.end();
	std::vector<ConfusionMatrixTuple> vcmt;
	
	// check the requested scores
	//Cans
	auto itE = std::find_if(start, end, [](char ch){return (ch == 'E' || ch == 'e');});
	auto itC = std::find_if(start, end, [](char ch){return (ch == 'C' || ch == 'c');});
	if(itC != end || itE != end)
	{
		if(itE != end)
		{
			findInteractionRef.set_isWobble_canonical(true);
			findInteractionQuery.set_isWobble_canonical(true);
		}
		vcmt.emplace_back(calcConfusionMatrixCansMatrix());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_Matrix_can, findInteractionQuery.m_Matrix_can));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
		
	//Wobbles
	auto itW = std::find_if(start, end, [](char ch){return (ch == 'W' || ch == 'w');});
	if(itW != end)
	{
		vcmt.emplace_back(calcConfusionMatrixWobblesMatrix());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_Matrix_wobble, findInteractionQuery.m_Matrix_wobble));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
	
	//NonCans
	auto itN = std::find_if(start, end, [](char ch){return (ch == 'N' || ch == 'n');});
	if(itN != end)
	{
		vcmt.emplace_back(calcConfusionMatrixNonCansMatrix());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_Matrix_noncan, findInteractionQuery.m_Matrix_noncan));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
	
	//All
	auto itA = std::find_if(start, end, [](char ch){return (ch == 'A' || ch == 'a');});
	if(itA != end)
	{
		vcmt.emplace_back(calcConfusionMatrixAllMatrix());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.m_Matrix_all, findInteractionQuery.m_Matrix_all));
	}
	else
	{
		vcmt.emplace_back(std::make_tuple(0,0,0,0));
		m_vMax_n_positives.push_back(0);
	}
	return vcmt;
}

std::map<std::string, double> SS::CompareStructures::calcSimilarityScores(ConfusionMatrixTuple cmt, int max_n_positives)
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
	simScoresMap.insert(std::make_pair("JBIndex", static_cast<double>(TP/(max_n_positives + 0.00005))));

	
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



void SS::CompareStructures::writeScores(std::vector<ConfusionMatrixTuple> vcmt, fs::path outputPath)
{
	std::ofstream outFile(outputPath.string().c_str());
	// writing the header to the file
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "F1," << "FM," << "J," 
		    << " PRECISION," << " RECALL," << " Specificity," << " BA," << " FOR,"
		    << "PT," << "CSI," << "MK," << "B" << std::endl;
	for(size_t i{0}; i < vcmt.size(); ++i)
	{
		auto [TP, TN, FP, FN] = vcmt[i];
		if(TP == 0 && TN == 0 && FP == 0 && FN == 0)
		{
			//outFile << vInteractionTypes1[i] << "-,-,-,-,-,-,-,-,-,-,-" << endl;
			continue;
		}
		std::map<std::string, double> simScoresMap = calcSimilarityScores(vcmt[i], m_vMax_n_positives[i]);
		auto itE = std::find_if(m_requestedInteractions.begin(), m_requestedInteractions.end(), [](char ch){return (ch == 'E' || ch == 'e');});
		if(itE == m_requestedInteractions.end())
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

void SS::CompareStructures::sepInteractions()
{
	findInteractionRef.mapSS2seq(findInteractionRef.m_SSMap);
	findInteractionQuery.mapSS2seq(findInteractionQuery.m_SSMap);
}

void SS::CompareStructures::readsequence(const fs::path& seqPath)
{
	findInteractionRef.init_sequence(seqPath);
	findInteractionQuery.init_sequence(seqPath);
	sepInteractions();
}


//#define _MAIN_
#ifdef _MAIN_
using namespace std;
int main()
{
	SS::CompareStructures cs;
	cs.m_requestedInteractions = "canwe";
	
	//cs.findInteractionRef.set_isWobble_canonical(false);
	//cs.findInteractionQuery.set_isWobble_canonical(false);
	
	//write the SSmaps of the ref and query
	cs.readStructures("samples/Cruciform.SS", "samples/Cruciform_local_all_minE-000001.ss_detected");
	cs.findInteractionRef.init_sequence("samples/Cruciform.fasta");
	cs.findInteractionQuery.init_sequence("samples/Cruciform.fasta");
	cs.sepInteractions();
	auto vcmt = cs.calcConfusionMatrixVector();
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
