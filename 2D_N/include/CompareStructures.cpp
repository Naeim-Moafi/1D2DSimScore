#include <algorithm> 
#include "CompareStructures.h"
#include <iomanip>
#include <fstream>

//constructors
CLARNA::CompareStructures::CompareStructures()
{
	requestedInteractions = "ABCNS";
}

CLARNA::CompareStructures::CompareStructures(std::string requestedInteractions)
: requestedInteractions(requestedInteractions) {}


// set_number_involved_faces_edges
void CLARNA::CompareStructures::set_number_involved_faces_edges(int number_involved_faces_edges)
{
	findInteractionRef.set_number_involved_faces_edges(number_involved_faces_edges);
	findInteractionQuery.set_number_involved_faces_edges(number_involved_faces_edges);
}

// set_isWobble_canonical
void CLARNA::CompareStructures::set_isWobble_canonical(bool isWobble_canonical)
{
	findInteractionRef.set_isWobble_canonical(isWobble_canonical);
	findInteractionQuery.set_isWobble_canonical(isWobble_canonical);
}

// readStructuresClaRNA
void CLARNA::CompareStructures::readStructures(const std::filesystem::path& refPath, const std::filesystem::path& queryPath, const std::filesystem::path& pdbPath)
{	
	auto vfct_pdb = findInteractionRef.readPDBFile(pdbPath);
	
	auto query = findInteractionQuery.readInputFile(queryPath);
	
	//Reference:
	auto vfct_clarna_ref = findInteractionRef.readInputFile(refPath);
	// all
	if(findInteractionRef.get_number_involved_faces_edges() == FIVE_EDGES)
	{
		findInteractionRef.addNonInteractedEdges(pdbPath, refPath);
	}
	// can
	findInteractionRef.sepCanBasePairs(vfct_clarna_ref, vfct_pdb);
	// noncan
	findInteractionRef.sepNonCanBasePairs(vfct_clarna_ref, vfct_pdb);
	// stacking
	findInteractionRef.sepStacking(vfct_clarna_ref, vfct_pdb);
	// wobbles
	findInteractionRef.sepWobble(vfct_clarna_ref, vfct_pdb);
	// all base pairs (canonical and noncanonical)
	findInteractionRef.sepBasePairs(vfct_clarna_ref, vfct_pdb);
	
	//Query:
	auto vfct_clarna_query = findInteractionQuery.readInputFile(queryPath);
	// all
	if(findInteractionQuery.get_number_involved_faces_edges() == FIVE_EDGES)
	{
		findInteractionQuery.addNonInteractedEdges(pdbPath, queryPath);
	}
	// can
	findInteractionQuery.sepCanBasePairs(vfct_clarna_query, vfct_pdb);
	// noncan
	findInteractionQuery.sepNonCanBasePairs(vfct_clarna_query, vfct_pdb);
	// stacking
	findInteractionQuery.sepStacking(vfct_clarna_query, vfct_pdb);
	// wobbles
	findInteractionQuery.sepWobble(vfct_clarna_query, vfct_pdb);
	// all base pairs (canonical and noncanonical)
	findInteractionQuery.sepBasePairs(vfct_clarna_query, vfct_pdb);
}

bool CLARNA::CompareStructures::isTP(FullClaRNATuple ftRef, FullClaRNATuple ftQuery)
{
	auto [_1R, _2R, _3R, _4R, secondChainIDRef, secondNuclNumberRef, secondNuclNameRef, secondEdgeRef, typeRef, weightRef] = ftRef;
	auto [_1Q, _2Q, _3Q, _4Q, secondChainIDQuery, secondNuclNumberQuery, secondNuclNameQuery, secondEdgeQuery, typeQuery, weightQuery] = ftQuery;
	// since the first nucleotide are the same 
	// if the second chainIDs
	// and the second nucleNumber
	// and the type of interaction (secondEdge and type) of reference and query were the same 
	// and also nucleNumber of the second nucleotide of the referene was not -1
	// it is true positive
	if(secondChainIDRef == secondChainIDQuery
	&& secondNuclNumberRef == secondNuclNumberQuery
	&& secondEdgeRef == secondEdgeQuery
	&& typeRef == typeQuery
	&& secondNuclNumberRef != -1)
	{
		return true;
	}
	
	return false;
}

// isTN
bool CLARNA::CompareStructures::isTN(FullClaRNATuple ftRef, FullClaRNATuple ftQuery)
{
	auto [_1R, _2R, _3R, _4R, secondChainIDRef, secondNuclNumberRef, secondNuclNameRef, secondEdgeRef, typeRef, weightRef] = ftRef;
	auto [_1Q, _2Q, _3Q, _4Q, secondChainIDQuery, secondNuclNumberQuery, secondNuclNameQuery, secondEdgeQuery, typeQuery, weightQuery] = ftQuery;		
	// if the second chainIDs
	// and the second nucleNumber
	// and the type of interaction of reference and query were the same
	// and also nucleNumber of the second nucleotide of the referene was -1
	// it is true negative
	if(secondChainIDRef == secondChainIDQuery
	&& secondNuclNumberRef == secondNuclNumberQuery
	&& secondNuclNumberRef == -1)
	{
		return true;
	}
	
	return false;
}

// isFN
bool CLARNA::CompareStructures::isFN(FullClaRNATuple ftRef, FullClaRNATuple ftQuery)
{
	auto [_1R, _2R, _3R, _4R, secondChainIDRef, secondNuclNumberRef, secondNuclNameRef, secondEdgeRef, typeRef, weightRef] = ftRef;
	auto [_1Q, _2Q, _3Q, _4Q, secondChainIDQuery, secondNuclNumberQuery, secondNuclNameQuery, secondEdgeQuery, typeQuery, weightQuery] = ftQuery;		
	// if the second chainIDs
	// or the second nucleNumbers
	// or the second edge
	// or the second type were not the same
	// and the second nucleNumber of the reference was not -1
	// it is false negative
	if((    secondChainIDRef != secondChainIDQuery
		 || secondNuclNumberRef != secondNuclNumberQuery
		 || secondEdgeRef != secondEdgeQuery
		 || typeRef != typeQuery)
	  )
	{
		if(secondNuclNumberRef != -1 && secondNuclNumberQuery == -1)
		{
			return true;
		}
	}
	
	return false;
}

// isFP
bool CLARNA::CompareStructures::isFP(FullClaRNATuple ftRef, FullClaRNATuple ftQuery)
{
	auto [_1R, _2R, _3R, _4R, secondChainIDRef, secondNuclNumberRef, secondNuclNameRef, secondEdgeRef, typeRef, weightRef] = ftRef;
	auto [_1Q, _2Q, _3Q, _4Q, secondChainIDQuery, secondNuclNumberQuery, secondNuclNameQuery, secondEdgeQuery, typeQuery, weightQuery] = ftQuery;		
	// if the second chainIDs
	// or the second nucleNumbers
	// or the second edge
	// or the second type were not the same
	// and the second nucleNumber of the reference was -1
	// it is false positive
	if((    secondChainIDRef != secondChainIDQuery
		 || secondNuclNumberRef != secondNuclNumberQuery
		 || secondEdgeRef != secondEdgeQuery
		 || typeRef != typeQuery)
	  )
	{
		if(secondNuclNumberRef == -1 || (secondNuclNumberRef != -1 && secondNuclNumberQuery != -1))
		{
			return true;
		}
	}
	
	return false;
}


int CLARNA::CompareStructures::find_max_n_positives(const std::vector<FullClaRNATuple>& ref, const std::vector<FullClaRNATuple>& query) const
{
	int np_ref = 0;
	int np_query = 0;
	
	for(auto fct : ref)
	{
		if(std::get<7>(fct) != "_")
		{
			++np_ref;
		}
	}
	
	
	for(auto fct : query)
	{
		if(std::get<7>(fct) != "_")
		{
			++np_query;
		}
	}
	
	if(np_query > np_ref)
	{
		return np_query;
	}
	
	return np_ref;
}

// calcConfusionMatrix
std::vector<ConfusionMatrixTuple> CLARNA::CompareStructures::calcConfusionMatrix(const std::filesystem::path& refPath, const std::filesystem::path& queryPath, const std::filesystem::path& pdbPath)
{
	//ConfusionMatrixTuple cmtAll, cmtCans, cmtNonCans, cmtStacks;
	std::vector<ConfusionMatrixTuple> vcmt;
	
	try
	{
		readStructures(refPath, queryPath, pdbPath);
	}
	catch(const std::invalid_argument& ex)
	{
		cerr << ex.what() << endl;
	}
	
	// checking A
	auto itA = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'A' or ch == 'a');});
	if(itA != cend(requestedInteractions))
	{
		if(findInteractionRef.get_number_involved_faces_edges() < FIVE_EDGES)
		{
			cout << "ERROR: Incorrect number of involved faces and edges for all interactions\n";
			cout << "you can use 5\n";
			exit(EXIT_FAILURE);
		}
		
		vcmt.push_back(calcConfusionMatrixForAll());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.vftAll, findInteractionQuery.vftAll));
	}
	else
	{
		vcmt.push_back(std::make_tuple(0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	
	// checking C
	auto itC = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'C' or ch == 'c');});
	auto itE = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' or ch == 'e');});
	if(itC != cend(requestedInteractions) || itE != cend(requestedInteractions))
	{	
		if(findInteractionRef.get_number_involved_faces_edges() == TWO_EDGES)
		{
			cout << "ERROR: Incorrect number of involved faces and edges for canonicals\n";
			cout << "you can use 1 or 3 or 5\n";
			exit(EXIT_FAILURE);
		}
		
		if(itE != cend(requestedInteractions))
		{	
			findInteractionRef.set_isWobble_canonical(true);
			findInteractionQuery.set_isWobble_canonical(true);
			cout << "Wobble considered as canonical\n";
		}
		vcmt.push_back(calcConfusionMatrixForCans());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.vftCans, findInteractionQuery.vftCans));
	}
	else
	{
		vcmt.push_back(std::make_tuple(0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	// checking N
	auto itN = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'N' or ch == 'n');});
	if(itN != cend(requestedInteractions))
	{	
		if(findInteractionRef.get_number_involved_faces_edges() == TWO_EDGES || findInteractionRef.get_number_involved_faces_edges() == ONE_EDGE)
		{
			cout << "ERROR: Incorrect number of involved faces and edges for noncanonical base pairs\n";
			cout << "you can use 3 or 5\n";
			exit(EXIT_FAILURE);
		}
		vcmt.push_back(calcConfusionMatrixForNonCans());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.vftNonCans, findInteractionQuery.vftNonCans));
	}
	else
	{
		vcmt.push_back(std::make_tuple(0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	
	// checking S
	auto itS = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'S' or ch == 's');});
	if(itS != cend(requestedInteractions))
	{
		if(findInteractionRef.get_number_involved_faces_edges() != TWO_EDGES && findInteractionRef.get_number_involved_faces_edges() != FIVE_EDGES)
		{
			cout << "ERROR: Incorrect number of involved faces and edges for stacking\n";
			cout << "you can use 2 or 5\n";
			exit(EXIT_FAILURE);
		}
		vcmt.push_back(calcConfusionMatrixForStacking());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.vftStacks, findInteractionQuery.vftStacks));
	}
	else
	{
		vcmt.push_back(std::make_tuple(0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	
	// checking W
	auto itW = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'W' or ch == 'w');});
	if(itW != cend(requestedInteractions))
	{
		if(findInteractionRef.get_number_involved_faces_edges() == TWO_EDGES)
		{
			cout << "ERROR: Incorrect number of involved faces and edges for wobbles\n";
			cout << "you can use 1 or 3 or 5\n";
			exit(EXIT_FAILURE);
		}
		vcmt.push_back(calcConfusionMatrixForWobbles());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.vftWobbles, findInteractionQuery.vftWobbles));
	}
	else
	{
		vcmt.push_back(std::make_tuple(0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	// checking B
	auto itB = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'B' or ch == 'b');});
	if(itB != cend(requestedInteractions))
	{
		if(findInteractionRef.get_number_involved_faces_edges() == TWO_EDGES || findInteractionRef.get_number_involved_faces_edges() == ONE_EDGE)
		{
			cout << "ERROR: Incorrect number of involved faces and edges for all base pairs\n";
			cout << "you can use 3 or 5\n";
			exit(EXIT_FAILURE);
		}
		vcmt.push_back(calcConfusionMatrixForBasePairs());
		m_vMax_n_positives.push_back(find_max_n_positives(findInteractionRef.vftBasePairs, findInteractionQuery.vftBasePairs));
	}
	else
	{
		vcmt.push_back(std::make_tuple(0, 0, 0, 0));
		m_vMax_n_positives.push_back(0);
	}
	
	return vcmt;
}


ConfusionMatrixTuple CLARNA::CompareStructures::calcConfusionMatrixForAll()
{	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < findInteractionRef.vftAll.size(); ++i)
	{
		if(isTP(findInteractionRef.vftAll[i], findInteractionQuery.vftAll[i]))
		{
			++TP;
		}
		
		if(isTN(findInteractionRef.vftAll[i], findInteractionQuery.vftAll[i]))
		{
			++TN;
		}
		
		if(isFN(findInteractionRef.vftAll[i], findInteractionQuery.vftAll[i]))
		{
			++FN;
		}
		
		if(isFP(findInteractionRef.vftAll[i], findInteractionQuery.vftAll[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN); 
}

// calcConfusionMatrixForCansClaRNA
ConfusionMatrixTuple CLARNA::CompareStructures::calcConfusionMatrixForCans()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < findInteractionRef.vftCans.size(); ++i)
	{
		if(isTP(findInteractionRef.vftCans[i], findInteractionQuery.vftCans[i]))
		{
			++TP;
		}
		
		if(isTN(findInteractionRef.vftCans[i], findInteractionQuery.vftCans[i]))
		{
			++TN;
		}
		
		if(isFN(findInteractionRef.vftCans[i], findInteractionQuery.vftCans[i]))
		{
			++FN;
		}
		
		if(isFP(findInteractionRef.vftCans[i], findInteractionQuery.vftCans[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

// calcConfusionMatrixForNonCans
ConfusionMatrixTuple CLARNA::CompareStructures::calcConfusionMatrixForNonCans()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < findInteractionRef.vftNonCans.size(); ++i)
	{
		if(isTP(findInteractionRef.vftNonCans[i], findInteractionQuery.vftNonCans[i]))
		{
			++TP;
		}
		
		if(isTN(findInteractionRef.vftNonCans[i], findInteractionQuery.vftNonCans[i]))
		{
			++TN;
		}
		
		if(isFN(findInteractionRef.vftNonCans[i], findInteractionQuery.vftNonCans[i]))
		{
			++FN;
		}
		
		if(isFP(findInteractionRef.vftNonCans[i], findInteractionQuery.vftNonCans[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

// calcConfusionMatrixForStacking
ConfusionMatrixTuple CLARNA::CompareStructures::calcConfusionMatrixForStacking()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < findInteractionRef.vftStacks.size(); ++i)
	{
		if(isTP(findInteractionRef.vftStacks[i], findInteractionQuery.vftStacks[i]))
		{
			++TP;
		}
		
		if(isTN(findInteractionRef.vftStacks[i], findInteractionQuery.vftStacks[i]))
		{
			++TN;
		}
		
		if(isFN(findInteractionRef.vftStacks[i], findInteractionQuery.vftStacks[i]))
		{
			++FN;
		}
		
		if(isFP(findInteractionRef.vftStacks[i], findInteractionQuery.vftStacks[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

// calcConfusionMatrixForWobbles
ConfusionMatrixTuple CLARNA::CompareStructures::calcConfusionMatrixForWobbles()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < findInteractionRef.vftBasePairs.size(); ++i)
	{
		if(isTP(findInteractionRef.vftWobbles[i], findInteractionQuery.vftWobbles[i]))
		{
			++TP;
		}
		
		if(isTN(findInteractionRef.vftWobbles[i], findInteractionQuery.vftWobbles[i]))
		{
			++TN;
		}
		
		if(isFN(findInteractionRef.vftWobbles[i], findInteractionQuery.vftWobbles[i]))
		{
			++FN;
		}
		
		if(isFP(findInteractionRef.vftWobbles[i], findInteractionQuery.vftWobbles[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}


// calcConfusionMatrixForBasePairs
ConfusionMatrixTuple CLARNA::CompareStructures::calcConfusionMatrixForBasePairs()
{
	
	int TP, TN, FP, FN;
	
	TP = TN = FP = FN = 0;
	
	// since the size of the structure must be the same
	// we can just use one loop for comparing structures
	for(size_t i = 0; i < findInteractionRef.vftBasePairs.size(); ++i)
	{
		if(isTP(findInteractionRef.vftBasePairs[i], findInteractionQuery.vftBasePairs[i]))
		{
			++TP;
		}
		
		if(isTN(findInteractionRef.vftBasePairs[i], findInteractionQuery.vftBasePairs[i]))
		{
			++TN;
		}
		
		if(isFN(findInteractionRef.vftBasePairs[i], findInteractionQuery.vftBasePairs[i]))
		{
			++FN;
		}
		
		if(isFP(findInteractionRef.vftBasePairs[i], findInteractionQuery.vftBasePairs[i]))
		{
			++FP;
		}
	}
	
	return std::make_tuple(TP, TN, FP, FN);
}

// calcSimilarityScores
// using nodiscard attribute to make sure the return value is not being ignored.
[[nodiscard]]std::map<std::string, double> CLARNA::CompareStructures::calcSimilarityScores(ConfusionMatrixTuple cmt, int max_n_positives)
{
	std::map<std::string, double> simScoresMap;
	// make a object of the SimilarityScores class with the variables of the confiusion matrix tuple
	SimilarityScores simScores(std::get<0>(cmt), std::get<1>(cmt), std::get<2>(cmt), std::get<3>(cmt));
	
	simScoresMap.insert(std::make_pair("MCC", simScores.calcMCC()));
	simScoresMap.insert(std::make_pair("FScore", simScores.calcFscore()));
	simScoresMap.insert(std::make_pair("FMIndex", simScores.calcFMIndex()));
	simScoresMap.insert(std::make_pair("JIndex", simScores.calcJIndex()));
	simScoresMap.insert(std::make_pair("Precision", simScores.calcPrecision()));
	simScoresMap.insert(std::make_pair("Recall", simScores.calcRecall()));
	simScoresMap.insert(std::make_pair("TP", std::get<0>(cmt)));
	simScoresMap.insert(std::make_pair("TN", std::get<1>(cmt)));
	simScoresMap.insert(std::make_pair("FP", std::get<2>(cmt)));
	simScoresMap.insert(std::make_pair("FN", std::get<3>(cmt)));
	simScoresMap.insert(std::make_pair("Specificity", simScores.calcSpecificty()));
	simScoresMap.insert(std::make_pair("BA", simScores.calcBA()));
	simScoresMap.insert(std::make_pair("FOR", simScores.calcFOR()));
	simScoresMap.insert(std::make_pair("PT", simScores.calcPT()));
	simScoresMap.insert(std::make_pair("CSI", simScores.calcCSI()));
	simScoresMap.insert(std::make_pair("MK", simScores.calcMK()));
	simScoresMap.insert(std::make_pair("JBIndex", static_cast<double>(std::get<0>(cmt)/(max_n_positives + 0.00005))));

	return simScoresMap;
}

void CLARNA::CompareStructures::writeScores(std::vector<ConfusionMatrixTuple> vcmtClaRNA, const std::filesystem::path& outputPath)
{
	using std::endl;
	std::ofstream outFile(outputPath.string());
	
	// writing the header to the file	
	outFile << "Interaction_Type," << "TP," << "TN," << "FP," << "FN," << "ALL," 
	    	<< " MCC,"  << "FSCORE," << "FM_INDEX," << "J_INDEX," 
		    << " PRECISION," << " RECALL," << " Specificity," << " BA," << " FOR,"
		    << "PT," << "CSI," << "MK," << "JBIndex" << std::endl;
	// cheking the requseted interactions
	// cheking for the calculation of the all;
	auto itA = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'A' or ch == 'a');});
	if(itA != cend(requestedInteractions))
	{
		if(findInteractionRef.get_number_involved_faces_edges() == FIVE_EDGES)
		{
			std::map<std::string, double> simScoresMapA = calcSimilarityScores(vcmtClaRNA[0], m_vMax_n_positives[0]);
			outFile << "All," << simScoresMapA["TP"] << "," << simScoresMapA["TN"] << "," << simScoresMapA["FP"] << "," 
					<< simScoresMapA["FN"] << "," << simScoresMapA["TP"] + simScoresMapA["TN"] + simScoresMapA["FP"] + simScoresMapA["FN"];
					if(simScoresMapA["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["MCC"];

					if(simScoresMapA["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["FScore"];

					if(simScoresMapA["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["FMIndex"];

					if(simScoresMapA["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["JIndex"];

					if(simScoresMapA["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["Precision"];

					if(simScoresMapA["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["Recall"];

					if(simScoresMapA["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["Specificity"];

					if(simScoresMapA["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["BA"];


					if(simScoresMapA["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["FOR"];

					if(simScoresMapA["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["PT"];

					if(simScoresMapA["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["CSI"];

					if(simScoresMapA["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["MK"];

					if(simScoresMapA["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapA["JBIndex"] << endl;
					
					
					//<< "," << simScoresMapA["MCC"] << "," 
					//<< simScoresMapA["FScore"] << "," << simScoresMapA["FMIndex"] << "," << simScoresMapA["JIndex"] << "," << simScoresMapA["Precision"] << ","	
					//<< simScoresMapA["Recall"] << "," << simScoresMapA["Specificity"] << "," << simScoresMapA["BA"] << "," << simScoresMapA["FOR"] << "," 
					//<< simScoresMapA["PT"] << "," << simScoresMapA["CSI"]  << "," << simScoresMapA["MK"] 
					//<< "," << simScoresMapA["JBIndex"] << endl;	
		}
		else
		{
			cout << "WARRNING: The number of the involved edges and faces are not enough for reporting the All types of interactions\n";
			outFile << "All, -, -, -, -, -, -, -, -, -, -, -\n";
		}
	}
	//else
	//{
		//outFile << "All, -, -, -, -, -, -, -, -, -, -, -\n";
	//}
	
	auto itC = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'C' or ch == 'c');});
	auto itE = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' or ch == 'e');});
	if(itC != cend(requestedInteractions) || itE != cend(requestedInteractions))
	{
		if(itE != cend(requestedInteractions))
		{
			std::map<std::string, double> simScoresMapC = calcSimilarityScores(vcmtClaRNA[1], m_vMax_n_positives[1]);
			outFile << "Cans + wobble," << simScoresMapC["TP"] << "," << simScoresMapC["TN"] << "," << simScoresMapC["FP"] << "," 
					<< simScoresMapC["FN"] << "," << simScoresMapC["TP"] + simScoresMapC["TN"] + simScoresMapC["FP"] + simScoresMapC["FN"];
					if(simScoresMapC["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["MCC"];

					if(simScoresMapC["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["FScore"];

					if(simScoresMapC["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["FMIndex"];

					if(simScoresMapC["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["JIndex"];

					if(simScoresMapC["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["Precision"];

					if(simScoresMapC["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["Recall"];

					if(simScoresMapC["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["Specificity"];

					if(simScoresMapC["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["BA"];


					if(simScoresMapC["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["FOR"];

					if(simScoresMapC["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["PT"];

					if(simScoresMapC["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["CSI"];

					if(simScoresMapC["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["MK"];

					if(simScoresMapC["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["JBIndex"] << endl;
					
					
					
					
					
					
					
					 //<< "," 
					//<< simScoresMapC["MCC"] << "," << simScoresMapC["FScore"] << "," 
					//<< simScoresMapC["FMIndex"] << "," << simScoresMapC["JIndex"] << "," << simScoresMapC["Precision"] << "," 
					//<< simScoresMapC["Recall"] << "," << simScoresMapC["Specificity"] << "," << simScoresMapC["BA"] << "," << simScoresMapC["FOR"] << "," 
					//<< simScoresMapC["PT"] << "," << simScoresMapC["CSI"]  << "," << simScoresMapC["MK"] 
				    //<< "," << simScoresMapC["JBIndex"] << endl;
			
		}
		else
		{
			std::map<std::string, double> simScoresMapC = calcSimilarityScores(vcmtClaRNA[1], m_vMax_n_positives[1]);
			outFile << "Cans," << simScoresMapC["TP"] << "," << simScoresMapC["TN"] << "," << simScoresMapC["FP"] << "," 
					<< simScoresMapC["FN"] << "," << simScoresMapC["TP"] + simScoresMapC["TN"] + simScoresMapC["FP"] + simScoresMapC["FN"];
					if(simScoresMapC["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["MCC"];

					if(simScoresMapC["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["FScore"];

					if(simScoresMapC["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["FMIndex"];

					if(simScoresMapC["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["JIndex"];

					if(simScoresMapC["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["Precision"];

					if(simScoresMapC["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["Recall"];

					if(simScoresMapC["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["Specificity"];

					if(simScoresMapC["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["BA"];


					if(simScoresMapC["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["FOR"];

					if(simScoresMapC["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["PT"];

					if(simScoresMapC["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["CSI"];

					if(simScoresMapC["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["MK"];

					if(simScoresMapC["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapC["JBIndex"] << endl;
					
					
					
					
					 //<< "," 
					//<< simScoresMapC["MCC"] << "," << simScoresMapC["FScore"] << "," 
					//<< simScoresMapC["FMIndex"] << "," << simScoresMapC["JIndex"] << "," << simScoresMapC["Precision"] << "," 
					//<< simScoresMapC["Recall"] << "," << simScoresMapC["Specificity"] << "," << simScoresMapC["BA"] << "," << simScoresMapC["FOR"] << "," 
					//<< simScoresMapC["PT"] << "," << simScoresMapC["CSI"]  << "," << simScoresMapC["MK"] 
				    //<< "," << simScoresMapC["JBIndex"] << endl;
		}
	}
	//else
	//{
		//outFile << "Cans, -, -, -, -, -, -, -, -, -, -, -\n";
	//}
	
	auto itN = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'N' or ch == 'n');});
	if(itN != cend(requestedInteractions))
	{
		std::map<std::string, double> simScoresMapN = calcSimilarityScores(vcmtClaRNA[2], m_vMax_n_positives[2]);
		outFile << "NonCans," << simScoresMapN["TP"] << "," << simScoresMapN["TN"] << "," << simScoresMapN["FP"] << ","
			    << simScoresMapN["FN"] << "," << simScoresMapN["TP"] + simScoresMapN["TN"] + simScoresMapN["FP"] + simScoresMapN["FN"];
					if(simScoresMapN["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["MCC"];

					if(simScoresMapN["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["FScore"];

					if(simScoresMapN["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["FMIndex"];

					if(simScoresMapN["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["JIndex"];

					if(simScoresMapN["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["Precision"];

					if(simScoresMapN["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["Recall"];

					if(simScoresMapN["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["Specificity"];

					if(simScoresMapN["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["BA"];


					if(simScoresMapN["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["FOR"];

					if(simScoresMapN["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["PT"];

					if(simScoresMapN["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["CSI"];

					if(simScoresMapN["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["MK"];

					if(simScoresMapN["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapN["JBIndex"] << endl;
			    
			    
			    
			    
			    
			     //<< "," 
			    //<< simScoresMapN["MCC"] << "," << simScoresMapN["FScore"] << "," 
				//<< simScoresMapN["FMIndex"] << "," << simScoresMapN["JIndex"] << "," << simScoresMapN["Precision"] << "," 
				//<< simScoresMapN["Recall"] << "," << simScoresMapN["Specificity"] << "," << simScoresMapN["BA"] << "," << simScoresMapN["FOR"] << "," 
				//<< simScoresMapN["PT"] << "," << simScoresMapN["CSI"]  << "," << simScoresMapN["MK"] 
				//<< "," << simScoresMapN["JBIndex"] << endl;	
	}
	//else
	//{
		//outFile << "NonCans, -, -, -, -, -, -, -, -, -, -, -\n";
	//}
	
	
	auto itW = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'W' or ch == 'w');});
	if(itW != cend(requestedInteractions))
	{
		std::map<std::string, double> simScoresMapW = calcSimilarityScores(vcmtClaRNA[4], m_vMax_n_positives[4]);
		outFile << "Wobbles," << simScoresMapW["TP"] << "," << simScoresMapW["TN"] << "," << simScoresMapW["FP"] << "," 
			    << simScoresMapW["FN"] << "," << simScoresMapW["TP"] + simScoresMapW["TN"] + simScoresMapW["FP"] + simScoresMapW["FN"];
					if(simScoresMapW["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["MCC"];

					if(simScoresMapW["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["FScore"];

					if(simScoresMapW["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["FMIndex"];

					if(simScoresMapW["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["JIndex"];

					if(simScoresMapW["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["Precision"];

					if(simScoresMapW["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["Recall"];

					if(simScoresMapW["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["Specificity"];

					if(simScoresMapW["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["BA"];


					if(simScoresMapW["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["FOR"];

					if(simScoresMapW["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["PT"];

					if(simScoresMapW["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["CSI"];

					if(simScoresMapW["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["MK"];

					if(simScoresMapW["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapW["JBIndex"] << endl;
			    
			    
			    
			    
			     //<< "," 
			    //<< simScoresMapW["MCC"] << "," << simScoresMapW["FScore"] << "," 
				//<< simScoresMapW["FMIndex"] << "," << simScoresMapW["JIndex"] << "," << simScoresMapW["Precision"] << "," 
				//<< simScoresMapW["Recall"] << "," << simScoresMapW["Specificity"] << "," << simScoresMapW["BA"] << "," << simScoresMapW["FOR"] << "," 
				//<< simScoresMapW["PT"] << "," << simScoresMapW["CSI"]  << "," << simScoresMapW["MK"] 
				//<< "," << simScoresMapW["JBIndex"] << endl;
	}
	//else
	//{
		//outFile << "Wobbles, -, -, -, -, -, -, -, -, -, -, -\n";
	//}
	
	auto itS = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'S' or ch == 's');});
	if(itS != cend(requestedInteractions))
	{
		std::map<std::string, double> simScoresMapS = calcSimilarityScores(vcmtClaRNA[3], m_vMax_n_positives[3]);
		outFile << "Stacks," << simScoresMapS["TP"] << "," << simScoresMapS["TN"] << "," << simScoresMapS["FP"] << "," 
			    << simScoresMapS["FN"] << "," << simScoresMapS["TP"] + simScoresMapS["TN"] + simScoresMapS["FP"] + simScoresMapS["FN"];
					if(simScoresMapS["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["MCC"];

					if(simScoresMapS["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["FScore"];

					if(simScoresMapS["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["FMIndex"];

					if(simScoresMapS["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["JIndex"];

					if(simScoresMapS["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["Precision"];

					if(simScoresMapS["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["Recall"];

					if(simScoresMapS["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["Specificity"];

					if(simScoresMapS["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["BA"];


					if(simScoresMapS["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["FOR"];

					if(simScoresMapS["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["PT"];

					if(simScoresMapS["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["CSI"];

					if(simScoresMapS["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["MK"];

					if(simScoresMapS["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapS["JBIndex"] << endl;
			    
			    
			    
			     //<< "," 
			    //<< simScoresMapS["MCC"] << "," << simScoresMapS["FScore"] << "," 
				//<< simScoresMapS["FMIndex"] << "," << simScoresMapS["JIndex"] << "," << simScoresMapS["Precision"] << "," 
				//<< simScoresMapS["Recall"] << "," << simScoresMapS["Specificity"] << "," << simScoresMapS["BA"] << "," << simScoresMapS["FOR"] << "," 
				//<< simScoresMapS["PT"] << "," << simScoresMapS["CSI"]  << "," << simScoresMapS["MK"] 
				//<< "," << simScoresMapS["JBIndex"] << endl;	
	}
	//else
	//{
		//outFile << "Stacks, -, -, -, -, -, -, -, -, -, -, -\n";
	//}
	
	
	auto itB = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'B' or ch == 'b');});
	if(itB != cend(requestedInteractions))
	{
		std::map<std::string, double> simScoresMapB = calcSimilarityScores(vcmtClaRNA[5], m_vMax_n_positives[5]);
		outFile << "BasePairs," << simScoresMapB["TP"] << "," << simScoresMapB["TN"] << "," << simScoresMapB["FP"] << "," 
			    << simScoresMapB["FN"] << "," << simScoresMapB["TP"] + simScoresMapB["TN"] + simScoresMapB["FP"] + simScoresMapB["FN"];
					if(simScoresMapB["MCC"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["MCC"];

					if(simScoresMapB["FScore"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["FScore"];

					if(simScoresMapB["FMIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["FMIndex"];

					if(simScoresMapB["JIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["JIndex"];

					if(simScoresMapB["Precision"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["Precision"];

					if(simScoresMapB["Recall"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["Recall"];

					if(simScoresMapB["Specificity"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["Specificity"];

					if(simScoresMapB["BA"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["BA"];


					if(simScoresMapB["FOR"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["FOR"];

					if(simScoresMapB["PT"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["PT"];

					if(simScoresMapB["CSI"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["CSI"];

					if(simScoresMapB["MK"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["MK"];

					if(simScoresMapB["JBIndex"] == -1.1) outFile << ",-";
					else outFile << "," << simScoresMapB["JBIndex"] << endl;
			    
			    
			    
			    
			    
			     //<< "," 
			    //<< simScoresMapB["MCC"] << "," << simScoresMapB["FScore"] << "," 
				//<< simScoresMapB["FMIndex"] << "," << simScoresMapB["JIndex"] << "," << simScoresMapB["Precision"] << "," 
				//<< simScoresMapB["Recall"] << "," << simScoresMapB["Specificity"] << "," << simScoresMapB["BA"] << "," << simScoresMapB["FOR"] << "," 
				//<< simScoresMapB["PT"] << "," << simScoresMapB["CSI"]  << "," << simScoresMapB["MK"] 
				//<< "," << simScoresMapB["JBIndex"] << endl;	
	}
	//else
	//{
		//outFile << "BasePairs, -, -, -, -, -, -, -, -, -, -, -\n";
	//}
	
	outFile.close();
}


//#define _MAIN_
#ifdef _MAIN_
using namespace std;
int main()
{
	CompareStructures cs("CAnsb");
	FindInteraction fi;
	
	cs.findInteractionRef.setWobble(true);
	cs.findInteractionQuery.setWobble(true);
	
	//ConfusionMatrixTuple cmtSS = cs.calcConfusionMatrixFromSS("smallAB1ss1.SS", "smallAB1ss2.SS");
	vector<ConfusionMatrixTuple> vcmtClaRNA = cs.calcConfusionMatrixFromClaRNA("ClaRNARef.out", "ClaRNAQuery.out", "ClaRNASeq.seq");
	for(const auto& cmt : vcmtClaRNA)
	{
		cout << "TP: " << get<0>(cmt) << " TN: " << get<1>(cmt) << " FP: " << get<2>(cmt) << " FN: " << get<3>(cmt) << endl;
	}
	
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
