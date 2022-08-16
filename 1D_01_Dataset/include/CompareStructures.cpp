#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

using std::endl;
const std::vector<std::string> vInteractionTypes1 = {"Cans,", "wobbles,", "NonCans,", "AllBasePairs,"};
const std::vector<std::string> vInteractionTypes2 = {"Cans + wobbles,", "wobbles,", "NonCans,", "AllBasePairs,"};


std::string toupper(const std::string& inputString)
{
	std::string s;
	for(auto c : inputString)
	{
		s += toupper(c);
	}
	return s;
}

std::vector<std::string> requestedScores_parser(const std::string& inputString)
{
	std::vector<std::string> v_scores;
	std::istringstream iss(toupper(inputString));
	std::string score;
	while(std::getline(iss, score, ','))
	{
		v_scores.emplace_back(score);
	}	
	
	return v_scores;
}

void Binary_all::CompareStructures::find_max_n_positives()
{
	
	for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
	{
		int n_positives = 0;
		for(size_t i_res { 0 }; i_res < findInteractions.m_vBinaryInterface[i].size(); ++i_res)
		{
			std::string res(1,findInteractions.m_vBinaryInterface[i][i_res]);
			if(res == "X" || res == "1")
			{
				++n_positives;
			}
		}
		if(n_positives > m_max_n_positives)
		{
			m_max_n_positives = n_positives;
		}
	}
}

void Binary_all::CompareStructures::init_binaries(const fs::path& path)
{
	findInteractions.init_binary(path);
}

ConfusionMatrixTuple Binary_all::CompareStructures::calcConfusionMatrix(const std::string& ref, const std::string& query)
{	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	for(size_t i_res { 0 }; i_res < ref.size(); ++i_res)
	{
		std::string ref_res(1,ref[i_res]);
		std::string query_res(1, query[i_res]);
		if(ref_res == query_res)
		{
			if(ref_res == "X" || ref_res == "1")
			{
				
				++TP;
			}
			else
			{
				++TN;
			}
		}
		else
		{
			if(ref_res == "X" || ref_res == "1")
			{
				++FN;
			}
			else
			{
				++FP;
			}
		}
	}
	
	
	return std::make_tuple(TP, TN, FP, FN);
}

ScoreMap Binary_all::CompareStructures::calcSimilarityScores()
{
ScoreMap scoreMap;
	ScoreMatrix sm;
	find_max_n_positives();
	auto v_scores = requestedScores_parser(m_requestedScores);
	cout << "Requested scores:\n";
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcMCC());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("MCC", sm));
	}
	
	//=============
	auto itFScore = std::find(v_scores.begin(), v_scores.end(), "FSCORE");
	if(itFScore != v_scores.end())
	{
		cout << "\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcFscore());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("FScore", sm));
	}
	
	//=============
	auto itJIndex = std::find(v_scores.begin(), v_scores.end(), "JINDEX");
	if(itJIndex != v_scores.end())
	{
		cout << "\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcJIndex());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JIndex", sm));
	}
	
	//=============
	auto itFMIndex = std::find(v_scores.begin(), v_scores.end(), "FMINDEX");
	if(itFMIndex != v_scores.end())
	{
		cout << "\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcFMIndex());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("FMIndex", sm));
	}
	
	//=============
	auto itRecall = std::find(v_scores.begin(), v_scores.end(), "RECALL");
	if(itRecall != v_scores.end())
	{
		cout << "\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcRecall());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("Recall", sm));
	}
	
	//=============
	auto itPrecision = std::find(v_scores.begin(), v_scores.end(), "PRECISION");
	if(itPrecision != v_scores.end())
	{
		cout << "\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcPrecision());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("Precision", sm));
	}
	
	//=============
	auto itSpecificity = std::find(v_scores.begin(), v_scores.end(), "SPECIFICITY");
	if(itSpecificity != v_scores.end())
	{
		cout << "\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcSpecificty());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("Specificity", sm));
	}
	
	//=============
	auto itBA = std::find(v_scores.begin(), v_scores.end(), "BA");
	if(itBA != v_scores.end())
	{
		cout << "\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcBA());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("BA", sm));
	}
	
	//=============
	auto itFOR = std::find(v_scores.begin(), v_scores.end(), "FOR");
	if(itFOR != v_scores.end())
	{
		cout << "\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcFOR());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("FOR", sm));
	}
	
	//=============
	auto itPT = std::find(v_scores.begin(), v_scores.end(), "PT");
	if(itPT != v_scores.end())
	{
		cout << "\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcPT());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("PT", sm));
	}
	
	//=============
	auto itCSI = std::find(v_scores.begin(), v_scores.end(), "CSI");
	if(itCSI != v_scores.end())
	{
		cout << "\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcCSI());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("CSI", sm));
	}
	
	//=============
	auto itMK = std::find(v_scores.begin(), v_scores.end(), "MK");
	if(itMK != v_scores.end())
	{
		cout << "\tMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(simScores.calcMK());
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("MK", sm));
	}	
		
	//=============
	auto itJBIndex = std::find(v_scores.begin(), v_scores.end(), "JBINDEX");
	if(itJBIndex != v_scores.end())
	{
		cout << "\tJBINDEX\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.m_vNames.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.m_vBinaryInterface[i], findInteractions.m_vBinaryInterface[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(m_max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}	
	
	
	return scoreMap;
}


void Binary_all::CompareStructures::writeScores(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScores();
	auto v_scores = requestedScores_parser(m_requestedScores);
	
	//=============	
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_MCC" << path.extension().string();
		std::ofstream outMCC(oss.str());
		outMCC << "#genesilico matrix .gsm\n";
		outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outMCC << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outMCC << " " << findInteractions.m_vNames[i];
		}
		outMCC << "\n";
		for(size_t i { 0 }; i < scoreMap["MCC"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["MCC"].size(); ++j)
			{
				outMCC << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["MCC"][i][j];
			}
			outMCC << "\n";
		}
		outMCC.close();
	}
	
	//=============
	auto itFScore = std::find(v_scores.begin(), v_scores.end(), "FSCORE");
	if(itFScore != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_FScore" << path.extension().string();
		std::ofstream outFScore(oss.str());
		outFScore << "#genesilico matrix .gsm\n";
		outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outFScore << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outFScore << " " << findInteractions.m_vNames[i];
		}
		outFScore << "\n";
		for(size_t i { 0 }; i < scoreMap["FScore"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["FScore"].size(); ++j)
			{
				outFScore << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["FScore"][i][j];
			}
			outFScore << "\n";
		}
		
		outFScore.close();
	}
		
	//=============
	auto itJIndex = std::find(v_scores.begin(), v_scores.end(), "JINDEX");
	if(itJIndex != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_JIndex" << path.extension().string();
		std::ofstream outJIndex(oss.str());
		outJIndex << "#genesilico matrix .gsm\n";
		outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outJIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outJIndex << " " << findInteractions.m_vNames[i];
		}
		outJIndex << "\n";
		for(size_t i { 0 }; i < scoreMap["JIndex"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["JIndex"].size(); ++j)
			{
				outJIndex << " "  << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["JIndex"][i][j];
			}
			outJIndex << "\n";
		}
		
		outJIndex.close();
	}
		
	//=============
	auto itFMIndex = std::find(v_scores.begin(), v_scores.end(), "FMINDEX");
	if(itFMIndex != v_scores.end())
	{
		std::ostringstream oss;
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_FMIndex" << path.extension().string();
		std::ofstream outFMIndex(oss.str());
		outFMIndex << "#genesilico matrix .gsm\n";
		outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outFMIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outFMIndex << " " << findInteractions.m_vNames[i];
		}
		outFMIndex << "\n";
		for(size_t i { 0 }; i < scoreMap["FMIndex"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["FMIndex"].size(); ++j)
			{
				outFMIndex << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["FMIndex"][i][j];
			}
			outFMIndex << "\n";
		}
		
		outFMIndex.close();
	}
		
	//=============
	auto itRecall = std::find(v_scores.begin(), v_scores.end(), "RECALL");
	if(itRecall != v_scores.end())
	{
		std::ostringstream oss;
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Recall" << path.extension().string();
		std::ofstream outRecall(oss.str());
		outRecall << "#genesilico matrix .gsm\n";
		outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outRecall << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outRecall << " " << findInteractions.m_vNames[i];
		}
		outRecall << "\n";
		for(size_t i { 0 }; i < scoreMap["Recall"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["Recall"].size(); ++j)
			{
				outRecall << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["Recall"][i][j];
			}
			outRecall << "\n";
		}
		
		outRecall.close();
	}
		
	//=============
	auto itPrecision = std::find(v_scores.begin(), v_scores.end(), "PRECISION");
	if(itPrecision != v_scores.end())
	{
		std::ostringstream oss;
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Precision" << path.extension().string();
		std::ofstream outPrecision(oss.str());
		outPrecision << "#genesilico matrix .gsm\n";
		outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outPrecision << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outPrecision << " " << findInteractions.m_vNames[i];
		}
		outPrecision << "\n";
		for(size_t i { 0 }; i < scoreMap["Precision"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["Precision"].size(); ++j)
			{
				outPrecision << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["Precision"][i][j];
			}
			outPrecision << "\n";
		}
		
		outPrecision.close();
	}
		
	//=============
	auto itSpecificity = std::find(v_scores.begin(), v_scores.end(), "SPECIFICITY");
	if(itSpecificity != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Specificity" << path.extension().string();
		std::ofstream outSpecificity(oss.str());
		outSpecificity << "#genesilico matrix .gsm\n";
		outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outSpecificity << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outSpecificity << " " << findInteractions.m_vNames[i];
		}
		outSpecificity << "\n";
		for(size_t i { 0 }; i < scoreMap["Specificity"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["Specificity"].size(); ++j)
			{
				outSpecificity << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["Specificity"][i][j];
			}
			outSpecificity << "\n";
		}
		
		outSpecificity.close();
	}
		
	//=============
	auto itBA = std::find(v_scores.begin(), v_scores.end(), "BA");
	if(itBA != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BA" << path.extension().string();
		std::ofstream outBA(oss.str());
		outBA << "#genesilico matrix .gsm\n";
		outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outBA << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outBA << " " << findInteractions.m_vNames[i];
		}
		outBA << "\n";
		for(size_t i { 0 }; i < scoreMap["BA"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["BA"].size(); ++j)
			{
				outBA << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["BA"][i][j];
			}
			outBA << "\n";
		}
		
		outBA.close();
	}
		
	//=============
	auto itFOR = std::find(v_scores.begin(), v_scores.end(), "FOR");
	if(itFOR != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_FOR" << path.extension().string();
		std::ofstream outFOR(oss.str());
		outFOR << "#genesilico matrix .gsm\n";
		outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outFOR << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outFOR << " " << findInteractions.m_vNames[i];
		}
		outFOR << "\n";
		for(size_t i { 0 }; i < scoreMap["FOR"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["FOR"].size(); ++j)
			{
				outFOR << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["FOR"][i][j];
			}
			outFOR << "\n";
		}
		
		outFOR.close();
	}
		
	//=============
	auto itPT = std::find(v_scores.begin(), v_scores.end(), "PT");
	if(itPT != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_PT" << path.extension().string();
		std::ofstream outPT(oss.str());
		outPT << "#genesilico matrix .gsm\n";
		outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outPT << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outPT << " " << findInteractions.m_vNames[i];
		}
		outPT << "\n";
		for(size_t i { 0 }; i < scoreMap["PT"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["PT"].size(); ++j)
			{
				outPT << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["PT"][i][j];
			}
			outPT << "\n";
		}
		
		outPT.close();
	}
		
	//=============
	auto itCSI = std::find(v_scores.begin(), v_scores.end(), "CSI");
	if(itCSI != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_CSI" << path.extension().string();
		std::ofstream outCSI(oss.str());
		outCSI << "#genesilico matrix .gsm\n";
		outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outCSI << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outCSI << " " << findInteractions.m_vNames[i];
		}
		outCSI << "\n";
		for(size_t i { 0 }; i < scoreMap["CSI"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["CSI"].size(); ++j)
			{
				outCSI << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["CSI"][i][j];
			}
			outCSI << "\n";
		}
		
		outCSI.close();
	}
		
	//=============
	auto itMK = std::find(v_scores.begin(), v_scores.end(), "MK");
	if(itMK != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_MK" << path.extension().string();
		std::ofstream outMK(oss.str());
		outMK << "#genesilico matrix .gsm\n";
		outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outMK << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outMK << " " << findInteractions.m_vNames[i];
		}
		outMK << "\n";
		for(size_t i { 0 }; i < scoreMap["MK"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["MK"].size(); ++j)
			{
				outMK << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["MK"][i][j];
			}
			outMK << "\n";
		}
		
		outMK.close();
	}
		
	//=============
	auto itJBIndex = std::find(v_scores.begin(), v_scores.end(), "JBINDEX");
	if(itJBIndex != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vNames.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vNames.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vNames[i];
		}
		outJBIndex << "\n";
		for(size_t i { 0 }; i < scoreMap["JBIndex"].size(); ++i)
		{
			for(size_t j { 0 }; j < scoreMap["JBIndex"].size(); ++j)
			{
				outJBIndex << " " << std::fixed << std::setw(5) << std::setprecision(3) << scoreMap["JBIndex"][i][j];
			}
			outJBIndex << "\n";
		}
		
		outJBIndex.close();
	}
	
	cout << "Your results are ready in " << path.parent_path() << " directory" << endl;
	
}

