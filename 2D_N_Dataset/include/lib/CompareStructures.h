/* 
 * Author:  Seyed Naeim Moafinejad
 */

/*
 * International Institute of Molecular and Cell Biology (IIMCB)
 * Copyright [2022] [IIMCB]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *        http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/


#pragma once

#include "FindInteractions.h"
#include "SimilarityScores.h"
#include <filesystem>
#include <map>

const int TRUE_POSITIVE = 1;
const int TRUE_NEGATIVE = 2;
const int FALSE_POSITIVE = 3;
const int FALSE_NEGATIVE = 4;

// define an alias type for a tuple which is represent the confusion matrix
using ConfusionMatrixTuple = std::tuple<int, int, int, int>; // (TP, TN, FP, FN);

namespace CLARNA{
	class CompareStructures
	{
		public:
			std::string requestedInteractions;
			FindInteractions findInteractionRef, findInteractionQuery;
			SimilarityScores simScores;
			
			// setter for some private parameter in the FindInteractions
			void set_number_involved_faces_edges(int number_involved_faces_edges);
			void set_isWobble_canonical(bool isWobble_canonical);

	
			// readStructures
			// by using the functions in the FindInteraction function read corresponding file
			// it gets paths of the reference and query structures.
			// and finally store the data in either vClaRNATuple;
			void readStructures(const std::filesystem::path& refPath, const std::filesystem::path& queryPath, const std::filesystem::path& pdbPath);
			
			//bool isTp(std::pair<int, int> pRef, std::pair<int, int> pQuery) const;
			bool isTP(FullClaRNATuple ftRef, FullClaRNATuple ftQuery);
			bool isTN(FullClaRNATuple ftRef, FullClaRNATuple ftQuery);
			bool isFP(FullClaRNATuple ftRef, FullClaRNATuple ftQuery);
			bool isFN(FullClaRNATuple ftRef, FullClaRNATuple ftQuery);
			
			// calcConfusionMatrix
			// gets the paths of the  reference and query structures
			// by using readStructures function will read structures and will compare the base-pairs			
			//			 if given base-pair with same edge and same type(cis and trans) was in both structure it would be 1 TP and 2 TN (for two other edges)
			//			 if there is no base pair in the same place in  both staructures at all it would be 3TN (for all edges)
			//			 if there is a base pair int the reference but not in query it would be 1 FN and 2 TN (for two other edges)
			//			 if there is not a given base-pair in reference but it was in query it would be 1 FP and 2 TN (for two other edges)
			[[nodiscard]] std::vector<ConfusionMatrixTuple> calcConfusionMatrix(const std::filesystem::path& refPath, const std::filesystem::path& queryPath, const std::filesystem::path& pdbPath);
			//// a helper function for above function
			//ConfusionMatrixTuple calcConfusionMatrixFromClaRNA(std::string requestedInteraction);
			
			// following function will calculate the confusion matrix
			// for the requested interactions
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrixForCans();
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrixForNonCans();
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrixForStacking();
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrixForBasePairs();
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrixForWobbles();
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrixForAll();
			
			// caclSimilarityScores gets a ConfusionMatrixTuple
			// and calculate all the score in the SimilarityScores class
			// and returns them as a map
			// in which the keys are the name of the scores and values are their value 
			std::map<std::string, double> calcSimilarityScores(ConfusionMatrixTuple cmt, int max_n_positives);
			
			// wirte the results into a file
			void writeScores(std::vector<ConfusionMatrixTuple> vcmtClaRNA, const std::filesystem::path& outputPath);
			
			CompareStructures();
			CompareStructures(std::string requestedIntercations);
			virtual ~CompareStructures() = default;
		private:
			std::vector<int> m_vMax_n_positives;
			int find_max_n_positives(const std::vector<FullClaRNATuple>& ref, const std::vector<FullClaRNATuple>& query) const;
	};
}
