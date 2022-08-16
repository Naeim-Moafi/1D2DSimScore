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
#include "lib/SimilarityScores.h"
#include "lib/CompareStructures.h"
#include <filesystem>
#include <map>

// define an alias type for a tuple which is represent the confusion matrix
//using ConfusionMatrixTuple = std::tuple<int, int, int, int>; // (TP, TN, FP, FN);

using ScoreMatrix = std::vector<std::vector<double>>;
using ScoreMap = std::map<std::string, ScoreMatrix>;

namespace CLARNA_all
{
	class CompareStructures : public CLARNA::CompareStructures
	{
		public:
			std::string m_requestedScores;
			FindInteractions findInteractions;
			std::string requestedInteractions;
			SimilarityScores simScores;
			
			// setter for some private parameter in the FindInteractions
			void set_number_involved_faces_edges(int number_involved_faces_edges);
			void set_isWobble_canonical(bool isWobble_canonical);
			
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrix(const std::vector<FullClaRNATuple>& vfct_ref, const std::vector<FullClaRNATuple>& vfct_query);
			
			// caclSimilarityScores gets a ConfusionMatrixTuple
			// and calculate all the score in the SimilarityScores class
			// and returns them as a map
			// in which the keys are the name of the scores and values are their value 
			[[nodiscard]] ScoreMap calcSimilarityScoresAll();
			[[nodiscard]] ScoreMap calcSimilarityScoresCans();
			[[nodiscard]] ScoreMap calcSimilarityScoresWobbles();
			[[nodiscard]] ScoreMap calcSimilarityScoresNonCans();
			[[nodiscard]] ScoreMap calcSimilarityScoresBasePairs();
			[[nodiscard]] ScoreMap calcSimilarityScoresStacks();
			
			
			// wirte the results into a file
			void writeScoresAll(const std::filesystem::path& path);
			void writeScoresCans(const std::filesystem::path& path);
			void writeScoresWobbles(const std::filesystem::path& path);
			void writeScoresNonCans(const std::filesystem::path& path);
			void writeScoresBasePairs(const std::filesystem::path& path);
			void writeScoresStacks(const std::filesystem::path& path);
			void writeScores(const std::filesystem::path& path);
			CompareStructures();
			CompareStructures(std::string requestedIntercations);
			virtual ~CompareStructures() = default;
		private:
			int find_max_n_positives(const std::vector<std::vector<FullClaRNATuple>>& v_vfct) const;
		
	};
}
