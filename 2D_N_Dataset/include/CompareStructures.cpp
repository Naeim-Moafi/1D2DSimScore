#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>


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



//constructors
CLARNA_all::CompareStructures::CompareStructures()
	: CLARNA::CompareStructures::CompareStructures() {}

CLARNA_all::CompareStructures::CompareStructures(std::string requestedInteractions)
: CLARNA::CompareStructures::CompareStructures(requestedInteractions) {}


// set_number_involved_faces_edges
void CLARNA_all::CompareStructures::set_number_involved_faces_edges(int number_involved_faces_edges)
{
	findInteractions.set_number_involved_faces_edges(number_involved_faces_edges);
}

// set_isWobble_canonical
void CLARNA_all::CompareStructures::set_isWobble_canonical(bool isWobble_canonical)
{
	findInteractions.set_isWobble_canonical(isWobble_canonical);
}


int CLARNA_all::CompareStructures::find_max_n_positives(const std::vector<std::vector<FullClaRNATuple>>& v_vfct) const
{
	int max_n_positives = 0;
	for(auto const& vfct : v_vfct)
	{
		int n_positives = 0;
		for(auto const& fct : vfct)
		{
			if(std::get<7>(fct) != "_")
			{
				++n_positives;
			}
		}
		if(max_n_positives < n_positives)
		{
			max_n_positives = n_positives;
		}
	}
	
	return max_n_positives;
}

ConfusionMatrixTuple CLARNA_all::CompareStructures::calcConfusionMatrix(const std::vector<FullClaRNATuple>& vfct_ref, const std::vector<FullClaRNATuple>& vfct_query)
{
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < vfct_ref.size(); ++i)
	{
		if(CLARNA::CompareStructures::isTP(vfct_ref[i], vfct_query[i]))
		{
			++TP;
		}
		
		if(CLARNA::CompareStructures::isTN(vfct_ref[i], vfct_query[i]))
		{
			++TN;
		}
		
		if(CLARNA::CompareStructures::isFN(vfct_ref[i], vfct_query[i]))
		{
			++FN;
		}
		
		if(CLARNA::CompareStructures::isFP(vfct_ref[i], vfct_query[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN); 
}

//calcSimilarityScores
ScoreMap CLARNA_all::CompareStructures::calcSimilarityScoresAll()
{
	ScoreMap scoreMap;
	ScoreMatrix sm;
	auto v_scores = requestedScores_parser(m_requestedScores);
	cout << "\tAll:\n";
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\t\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
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
		cout << "\t\tJBIndex\n";
		int max_n_positives = find_max_n_positives(findInteractions.v_vftAll);
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftAll.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftAll.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftAll[i], findInteractions.v_vftAll[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}	
	
	
	return scoreMap;
}


ScoreMap CLARNA_all::CompareStructures::calcSimilarityScoresCans()
{
	ScoreMap scoreMap;
	ScoreMatrix sm;
	auto v_scores = requestedScores_parser(m_requestedScores);
		
	auto itE = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' or ch == 'e');});
	if(itE != cend(requestedInteractions)) 
	{
		set_isWobble_canonical(true);
		cout << "\tCans+wobble:\n";
	}
	else
	{
		cout << "\tCans:\n";
	}
		
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\t\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
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
		cout << "\t\tJBIndex\n";
		int max_n_positives = find_max_n_positives(findInteractions.v_vftCans);
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftCans[i], findInteractions.v_vftCans[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}
	
	
	return scoreMap;
}


ScoreMap CLARNA_all::CompareStructures::calcSimilarityScoresWobbles()
{
	ScoreMap scoreMap;
	ScoreMatrix sm;
	auto v_scores = requestedScores_parser(m_requestedScores);
	cout << "\tWobbles:\n";
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\t\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\ttMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
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
		cout << "\t\tJBIndex\n";
		int max_n_positives = find_max_n_positives(findInteractions.v_vftWobbles);
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftWobbles.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftWobbles.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftWobbles[i], findInteractions.v_vftWobbles[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}
	
	
	return scoreMap;
}


ScoreMap CLARNA_all::CompareStructures::calcSimilarityScoresNonCans()
{
	ScoreMap scoreMap;
	ScoreMatrix sm;
	auto v_scores = requestedScores_parser(m_requestedScores);
	cout << "\tNonCans:\n";
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\t\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\ttMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
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
		cout << "\t\tJBIndex\n";
		int max_n_positives = find_max_n_positives(findInteractions.v_vftNonCans);
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftNonCans.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftNonCans.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftNonCans[i], findInteractions.v_vftNonCans[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}
	
	
	return scoreMap;
}


ScoreMap CLARNA_all::CompareStructures::calcSimilarityScoresBasePairs()
{
	ScoreMap scoreMap;
	ScoreMatrix sm;
	auto v_scores = requestedScores_parser(m_requestedScores);
	cout << "\tBasePairs:\n";
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\t\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
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
		cout << "\t\tJBIndex\n";
		int max_n_positives = find_max_n_positives(findInteractions.v_vftBasePairs);
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftBasePairs.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftBasePairs.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftBasePairs[i], findInteractions.v_vftBasePairs[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}	
	
	
	return scoreMap;
}


ScoreMap CLARNA_all::CompareStructures::calcSimilarityScoresStacks()
{
	ScoreMap scoreMap;
	ScoreMatrix sm;
	auto v_scores = requestedScores_parser(m_requestedScores);
	cout << "\tStacks:\n";
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		cout << "\t\tMCC\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tFScore\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tJIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tFMIndex\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tRecall\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tPrecision\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tSpecificity\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tBA\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tFOR\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tPT\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tCSI\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tMK\n";
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
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
		cout << "\t\tJBIndex\n";
		int max_n_positives = find_max_n_positives(findInteractions.v_vftStacks);
		sm.clear();
		for(size_t i { 0 }; i < findInteractions.v_vftStacks.size(); ++i)
		{
			std::vector<double> tmp_vec;
			for(size_t j { 0 }; j < findInteractions.v_vftStacks.size(); ++j)
			{
				auto [TP, TN, FP, FN] = calcConfusionMatrix(findInteractions.v_vftStacks[i], findInteractions.v_vftStacks[j]);
				simScores.setCategories(TP, TN, FP, FN);
				tmp_vec.push_back(static_cast<double>(TP/(max_n_positives + 0.00005)));
			}
			sm.emplace_back(tmp_vec);
		}
		scoreMap.insert(std::make_pair("JBIndex", sm));
	}	
	
	
	return scoreMap;
}

void CLARNA_all::CompareStructures::writeScoresAll(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScoresAll();
	auto v_scores = requestedScores_parser(m_requestedScores);
	//=============	
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_MCC" << path.extension().string();
		std::ofstream outMCC(oss.str());
		outMCC << "#genesilico matrix .gsm\n";
		outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMCC << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMCC << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_FScore" << path.extension().string();
		std::ofstream outFScore(oss.str());
		outFScore << "#genesilico matrix .gsm\n";
		outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFScore << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFScore << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_JIndex" << path.extension().string();
		std::ofstream outJIndex(oss.str());
		outJIndex << "#genesilico matrix .gsm\n";
		outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_All_FMIndex" << path.extension().string();
		std::ofstream outFMIndex(oss.str());
		outFMIndex << "#genesilico matrix .gsm\n";
		outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFMIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_All_Recall" << path.extension().string();
		std::ofstream outRecall(oss.str());
		outRecall << "#genesilico matrix .gsm\n";
		outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outRecall << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outRecall << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_All_Precision" << path.extension().string();
		std::ofstream outPrecision(oss.str());
		outPrecision << "#genesilico matrix .gsm\n";
		outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPrecision << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPrecision << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_Specificity" << path.extension().string();
		std::ofstream outSpecificity(oss.str());
		outSpecificity << "#genesilico matrix .gsm\n";
		outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outSpecificity << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_BA" << path.extension().string();
		std::ofstream outBA(oss.str());
		outBA << "#genesilico matrix .gsm\n";
		outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outBA << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outBA << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_FOR" << path.extension().string();
		std::ofstream outFOR(oss.str());
		outFOR << "#genesilico matrix .gsm\n";
		outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFOR << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFOR << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_PT" << path.extension().string();
		std::ofstream outPT(oss.str());
		outPT << "#genesilico matrix .gsm\n";
		outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPT << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPT << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_CSI" << path.extension().string();
		std::ofstream outCSI(oss.str());
		outCSI << "#genesilico matrix .gsm\n";
		outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outCSI << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outCSI << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_MK" << path.extension().string();
		std::ofstream outMK(oss.str());
		outMK << "#genesilico matrix .gsm\n";
		outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMK << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMK << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_All_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
}

void CLARNA_all::CompareStructures::writeScoresCans(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScoresCans();
	auto v_scores = requestedScores_parser(m_requestedScores);	
	auto itE = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' or ch == 'e');});
	if(itE != cend(requestedInteractions)) 
	{
		//=============	
		auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
		if(itMCC != v_scores.end())
		{
			std::ostringstream oss;
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_MCC" << path.extension().string();
			std::ofstream outMCC(oss.str());
			outMCC << "#genesilico matrix .gsm\n";
			outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outMCC << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outMCC << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_FScore" << path.extension().string();
			std::ofstream outFScore(oss.str());
			outFScore << "#genesilico matrix .gsm\n";
			outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outFScore << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outFScore << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_JIndex" << path.extension().string();
			std::ofstream outJIndex(oss.str());
			outJIndex << "#genesilico matrix .gsm\n";
			outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outJIndex << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outJIndex << " " << findInteractions.m_vStructures_name[i];
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
			oss  << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_FMIndex" << path.extension().string();
			std::ofstream outFMIndex(oss.str());
			outFMIndex << "#genesilico matrix .gsm\n";
			outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outFMIndex << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
			oss  << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_Recall" << path.extension().string();
			std::ofstream outRecall(oss.str());
			outRecall << "#genesilico matrix .gsm\n";
			outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outRecall << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outRecall << " " << findInteractions.m_vStructures_name[i];
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
			oss  << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_Precision" << path.extension().string();
			std::ofstream outPrecision(oss.str());
			outPrecision << "#genesilico matrix .gsm\n";
			outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outPrecision << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outPrecision << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_Specificity" << path.extension().string();
			std::ofstream outSpecificity(oss.str());
			outSpecificity << "#genesilico matrix .gsm\n";
			outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outSpecificity << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_BA" << path.extension().string();
			std::ofstream outBA(oss.str());
			outBA << "#genesilico matrix .gsm\n";
			outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outBA << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outBA << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_FOR" << path.extension().string();
			std::ofstream outFOR(oss.str());
			outFOR << "#genesilico matrix .gsm\n";
			outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outFOR << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outFOR << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_PT" << path.extension().string();
			std::ofstream outPT(oss.str());
			outPT << "#genesilico matrix .gsm\n";
			outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outPT << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outPT << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_CSI" << path.extension().string();
			std::ofstream outCSI(oss.str());
			outCSI << "#genesilico matrix .gsm\n";
			outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outCSI << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outCSI << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_MK" << path.extension().string();
			std::ofstream outMK(oss.str());
			outMK << "#genesilico matrix .gsm\n";
			outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outMK << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outMK << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans+wobble_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
	}
	else
	{
		//=============	
		auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
		if(itMCC != v_scores.end())
		{
			std::ostringstream oss;
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans_MCC" << path.extension().string();
			std::ofstream outMCC(oss.str());
			outMCC << "#genesilico matrix .gsm\n";
			outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outMCC << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outMCC << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans_FScore" << path.extension().string();
			std::ofstream outFScore(oss.str());
			outFScore << "#genesilico matrix .gsm\n";
			outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outFScore << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outFScore << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans_JIndex" << path.extension().string();
			std::ofstream outJIndex(oss.str());
			outJIndex << "#genesilico matrix .gsm\n";
			outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outJIndex << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outJIndex << " " << findInteractions.m_vStructures_name[i];
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
			oss  << path.parent_path().string() << "/" << path.stem().string() << "_Cans_FMIndex" << path.extension().string();
			std::ofstream outFMIndex(oss.str());
			outFMIndex << "#genesilico matrix .gsm\n";
			outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outFMIndex << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
			oss  << path.parent_path().string() << "/" << path.stem().string() << "_Cans__Recall" << path.extension().string();
			std::ofstream outRecall(oss.str());
			outRecall << "#genesilico matrix .gsm\n";
			outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outRecall << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outRecall << " " << findInteractions.m_vStructures_name[i];
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
			oss  << path.parent_path().string() << "/" << path.stem().string() << "_Cans__Precision" << path.extension().string();
			std::ofstream outPrecision(oss.str());
			outPrecision << "#genesilico matrix .gsm\n";
			outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outPrecision << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outPrecision << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans__Specificity" << path.extension().string();
			std::ofstream outSpecificity(oss.str());
			outSpecificity << "#genesilico matrix .gsm\n";
			outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outSpecificity << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans__BA" << path.extension().string();
			std::ofstream outBA(oss.str());
			outBA << "#genesilico matrix .gsm\n";
			outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outBA << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outBA << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans__FOR" << path.extension().string();
			std::ofstream outFOR(oss.str());
			outFOR << "#genesilico matrix .gsm\n";
			outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outFOR << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outFOR << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans__PT" << path.extension().string();
			std::ofstream outPT(oss.str());
			outPT << "#genesilico matrix .gsm\n";
			outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outPT << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outPT << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans__CSI" << path.extension().string();
			std::ofstream outCSI(oss.str());
			outCSI << "#genesilico matrix .gsm\n";
			outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outCSI << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outCSI << " " << findInteractions.m_vStructures_name[i];
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
			oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans__MK" << path.extension().string();
			std::ofstream outMK(oss.str());
			outMK << "#genesilico matrix .gsm\n";
			outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
			outMK << "#Names:";
			for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
			{
				outMK << " " << findInteractions.m_vStructures_name[i];
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
	}
		
	//=============
	auto itJBIndex = std::find(v_scores.begin(), v_scores.end(), "JBINDEX");
	if(itJBIndex != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Cans_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
	
}



void CLARNA_all::CompareStructures::writeScoresWobbles(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScoresWobbles();
	auto v_scores = requestedScores_parser(m_requestedScores);
	
	//=============	
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_MCC" << path.extension().string();
		std::ofstream outMCC(oss.str());
		outMCC << "#genesilico matrix .gsm\n";
		outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMCC << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMCC << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_FScore" << path.extension().string();
		std::ofstream outFScore(oss.str());
		outFScore << "#genesilico matrix .gsm\n";
		outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFScore << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFScore << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_JIndex" << path.extension().string();
		std::ofstream outJIndex(oss.str());
		outJIndex << "#genesilico matrix .gsm\n";
		outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_FMIndex" << path.extension().string();
		std::ofstream outFMIndex(oss.str());
		outFMIndex << "#genesilico matrix .gsm\n";
		outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFMIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_Recall" << path.extension().string();
		std::ofstream outRecall(oss.str());
		outRecall << "#genesilico matrix .gsm\n";
		outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outRecall << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outRecall << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_Precision" << path.extension().string();
		std::ofstream outPrecision(oss.str());
		outPrecision << "#genesilico matrix .gsm\n";
		outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPrecision << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPrecision << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_Specificity" << path.extension().string();
		std::ofstream outSpecificity(oss.str());
		outSpecificity << "#genesilico matrix .gsm\n";
		outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outSpecificity << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_BA" << path.extension().string();
		std::ofstream outBA(oss.str());
		outBA << "#genesilico matrix .gsm\n";
		outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outBA << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outBA << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_FOR" << path.extension().string();
		std::ofstream outFOR(oss.str());
		outFOR << "#genesilico matrix .gsm\n";
		outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFOR << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFOR << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_PT" << path.extension().string();
		std::ofstream outPT(oss.str());
		outPT << "#genesilico matrix .gsm\n";
		outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPT << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPT << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_CSI" << path.extension().string();
		std::ofstream outCSI(oss.str());
		outCSI << "#genesilico matrix .gsm\n";
		outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outCSI << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outCSI << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_MK" << path.extension().string();
		std::ofstream outMK(oss.str());
		outMK << "#genesilico matrix .gsm\n";
		outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMK << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMK << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Wobbles_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
}



void CLARNA_all::CompareStructures::writeScoresNonCans(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScoresNonCans();
	auto v_scores = requestedScores_parser(m_requestedScores);
	
	//=============	
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_MCC" << path.extension().string();
		std::ofstream outMCC(oss.str());
		outMCC << "#genesilico matrix .gsm\n";
		outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMCC << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMCC << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_FScore" << path.extension().string();
		std::ofstream outFScore(oss.str());
		outFScore << "#genesilico matrix .gsm\n";
		outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFScore << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFScore << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_JIndex" << path.extension().string();
		std::ofstream outJIndex(oss.str());
		outJIndex << "#genesilico matrix .gsm\n";
		outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_FMIndex" << path.extension().string();
		std::ofstream outFMIndex(oss.str());
		outFMIndex << "#genesilico matrix .gsm\n";
		outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFMIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_Recall" << path.extension().string();
		std::ofstream outRecall(oss.str());
		outRecall << "#genesilico matrix .gsm\n";
		outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outRecall << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outRecall << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_Precision" << path.extension().string();
		std::ofstream outPrecision(oss.str());
		outPrecision << "#genesilico matrix .gsm\n";
		outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPrecision << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPrecision << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_Specificity" << path.extension().string();
		std::ofstream outSpecificity(oss.str());
		outSpecificity << "#genesilico matrix .gsm\n";
		outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outSpecificity << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_BA" << path.extension().string();
		std::ofstream outBA(oss.str());
		outBA << "#genesilico matrix .gsm\n";
		outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outBA << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outBA << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_FOR" << path.extension().string();
		std::ofstream outFOR(oss.str());
		outFOR << "#genesilico matrix .gsm\n";
		outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFOR << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFOR << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_PT" << path.extension().string();
		std::ofstream outPT(oss.str());
		outPT << "#genesilico matrix .gsm\n";
		outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPT << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPT << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_CSI" << path.extension().string();
		std::ofstream outCSI(oss.str());
		outCSI << "#genesilico matrix .gsm\n";
		outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outCSI << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outCSI << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans__MK" << path.extension().string();
		std::ofstream outMK(oss.str());
		outMK << "#genesilico matrix .gsm\n";
		outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMK << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMK << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_NonCans_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
}


void CLARNA_all::CompareStructures::writeScoresBasePairs(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScoresBasePairs();
	auto v_scores = requestedScores_parser(m_requestedScores);
	
	//=============	
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_MCC" << path.extension().string();
		std::ofstream outMCC(oss.str());
		outMCC << "#genesilico matrix .gsm\n";
		outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMCC << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMCC << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_FScore" << path.extension().string();
		std::ofstream outFScore(oss.str());
		outFScore << "#genesilico matrix .gsm\n";
		outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFScore << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFScore << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_JIndex" << path.extension().string();
		std::ofstream outJIndex(oss.str());
		outJIndex << "#genesilico matrix .gsm\n";
		outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_FMIndex" << path.extension().string();
		std::ofstream outFMIndex(oss.str());
		outFMIndex << "#genesilico matrix .gsm\n";
		outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFMIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_Recall" << path.extension().string();
		std::ofstream outRecall(oss.str());
		outRecall << "#genesilico matrix .gsm\n";
		outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outRecall << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outRecall << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_Precision" << path.extension().string();
		std::ofstream outPrecision(oss.str());
		outPrecision << "#genesilico matrix .gsm\n";
		outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPrecision << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPrecision << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_Specificity" << path.extension().string();
		std::ofstream outSpecificity(oss.str());
		outSpecificity << "#genesilico matrix .gsm\n";
		outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outSpecificity << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_BA" << path.extension().string();
		std::ofstream outBA(oss.str());
		outBA << "#genesilico matrix .gsm\n";
		outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outBA << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outBA << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_FOR" << path.extension().string();
		std::ofstream outFOR(oss.str());
		outFOR << "#genesilico matrix .gsm\n";
		outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFOR << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFOR << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_PT" << path.extension().string();
		std::ofstream outPT(oss.str());
		outPT << "#genesilico matrix .gsm\n";
		outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPT << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPT << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_CSI" << path.extension().string();
		std::ofstream outCSI(oss.str());
		outCSI << "#genesilico matrix .gsm\n";
		outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outCSI << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outCSI << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_MK" << path.extension().string();
		std::ofstream outMK(oss.str());
		outMK << "#genesilico matrix .gsm\n";
		outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMK << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMK << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_BasePairs_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
}



void CLARNA_all::CompareStructures::writeScoresStacks(const std::filesystem::path& path)
{
	auto scoreMap = calcSimilarityScoresStacks();
	auto v_scores = requestedScores_parser(m_requestedScores);
	
	//=============	
	auto itMCC = std::find(v_scores.begin(), v_scores.end(), "MCC");
	if(itMCC != v_scores.end())
	{
		std::ostringstream oss;
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_MCC" << path.extension().string();
		std::ofstream outMCC(oss.str());
		outMCC << "#genesilico matrix .gsm\n";
		outMCC << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMCC << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMCC << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_FScore" << path.extension().string();
		std::ofstream outFScore(oss.str());
		outFScore << "#genesilico matrix .gsm\n";
		outFScore << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFScore << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFScore << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_JIndex" << path.extension().string();
		std::ofstream outJIndex(oss.str());
		outJIndex << "#genesilico matrix .gsm\n";
		outJIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_FMIndex" << path.extension().string();
		std::ofstream outFMIndex(oss.str());
		outFMIndex << "#genesilico matrix .gsm\n";
		outFMIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFMIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFMIndex << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_Recall" << path.extension().string();
		std::ofstream outRecall(oss.str());
		outRecall << "#genesilico matrix .gsm\n";
		outRecall << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outRecall << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outRecall << " " << findInteractions.m_vStructures_name[i];
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
		oss  << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_Precision" << path.extension().string();
		std::ofstream outPrecision(oss.str());
		outPrecision << "#genesilico matrix .gsm\n";
		outPrecision << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPrecision << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPrecision << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_Specificity" << path.extension().string();
		std::ofstream outSpecificity(oss.str());
		outSpecificity << "#genesilico matrix .gsm\n";
		outSpecificity << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outSpecificity << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outSpecificity << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_BA" << path.extension().string();
		std::ofstream outBA(oss.str());
		outBA << "#genesilico matrix .gsm\n";
		outBA << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outBA << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outBA << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_FOR" << path.extension().string();
		std::ofstream outFOR(oss.str());
		outFOR << "#genesilico matrix .gsm\n";
		outFOR << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outFOR << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outFOR << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_PT" << path.extension().string();
		std::ofstream outPT(oss.str());
		outPT << "#genesilico matrix .gsm\n";
		outPT << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outPT << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outPT << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_CSI" << path.extension().string();
		std::ofstream outCSI(oss.str());
		outCSI << "#genesilico matrix .gsm\n";
		outCSI << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outCSI << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outCSI << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_MK" << path.extension().string();
		std::ofstream outMK(oss.str());
		outMK << "#genesilico matrix .gsm\n";
		outMK << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outMK << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outMK << " " << findInteractions.m_vStructures_name[i];
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
		oss << path.parent_path().string() << "/" << path.stem().string() << "_Stacks_JBIndex" << path.extension().string();
		std::ofstream outJBIndex(oss.str());
		outJBIndex << "#genesilico matrix .gsm\n";
		outJBIndex << "#Format	DataType=Similarity DataFormat=Square Number=" << findInteractions.m_vStructures_name.size() << "\n";
		outJBIndex << "#Names:";
		for(size_t i { 0 }; i < findInteractions.m_vStructures_name.size(); ++i)
		{
			outJBIndex << " " << findInteractions.m_vStructures_name[i];
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
}


void CLARNA_all::CompareStructures::writeScores(const std::filesystem::path& path)
{
	auto itA = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'A' or ch == 'a');});
	if(itA != cend(requestedInteractions)) 
	{
		writeScoresAll(path);
	}
	
	
	auto itC = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'C' or ch == 'c');});
	auto itE = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' or ch == 'e');});
	if(itC != cend(requestedInteractions) || itE != cend(requestedInteractions)) 
	{
		writeScoresCans(path);
	}
	
	
	 auto itW = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'W' or ch == 'w');});
	if(itW != cend(requestedInteractions)) 
	{
		writeScoresWobbles(path);
	}
	
	
	auto itN = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'N' or ch == 'n');});
	if(itN != cend(requestedInteractions)) 
	{
		writeScoresNonCans(path);
	}
	
	
	auto itB = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'B' or ch == 'b');});
	if(itB != cend(requestedInteractions)) 
	{
		writeScoresBasePairs(path);
	}
	
	
	auto itS = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'S' or ch == 's');});
	if(itS != cend(requestedInteractions)) 
	{
		writeScoresStacks(path);
	}
}
