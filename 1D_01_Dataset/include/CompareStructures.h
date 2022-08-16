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
#include <map>
#include "lib/SimilarityScores.h"
#include <utility>

const int TRUE_POSITIVE = 1;
const int TRUE_NEGATIVE = 2;
const int FALSE_POSITIVE = 3;
const int FALSE_NEGATIVE = 4;

// define an alias type for a tuple which is represent the confusion matrix
using ConfusionMatrixTuple = std::tuple<int, int, int, int>; // (TP, TN, FP, FN);
namespace fs = std::filesystem;
using std::cout;
using std::cerr;
using std::endl;

using ScoreMatrix = std::vector<std::vector<double>>;
using ScoreMap = std::map<std::string, ScoreMatrix>;

namespace Binary_all
{
	class CompareStructures
	{
		public:
			std::string m_requestedScores;
			FindInteractions findInteractions;
			SimilarityScores simScores;
			ScoreMatrix scoreMatrix;
			
			// readStructures(
			// by using the functions in the FindInteraction class read corresponding file
			// it gets paths of the reference and query structures.
			// and finally store the data in either vClaRNATuple or vSS;
			void init_binaries(const fs::path& Path);
			
			
			
			// calcConfusionMatrix
			// gets two string for reference and query binary interaface interaction
			// and  calculate the confusion matrix's category.
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrix(const std::string& ref, const std::string& query);
			
			// caclSimilarityScores gets a ConfusionMatrixTuple
			// and calculate all the score in the SimilarityScores class
			// and returns them as a map
			// in which the keys are the name of the scores and values are their value 
			[[nodiscard]] ScoreMap calcSimilarityScores();
			
			
			// wirte the results into a file
			void writeScores(const std::filesystem::path& path);
			
			CompareStructures() = default;
			virtual ~CompareStructures() = default;
			
			private:
				int m_max_n_positives = 0;
				void find_max_n_positives();
		
	};
}
