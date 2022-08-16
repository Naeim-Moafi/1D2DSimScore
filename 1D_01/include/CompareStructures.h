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


namespace Binary
{
	class CompareStructures
	{
		public:
			FindInteractions findInteractionRef, findInteractionQuery;
			SimilarityScores simScores;
			int m_max_n_positives = 0;
			
			// readStructures(
			// by using the functions in the FindInteraction class read corresponding file
			// it gets paths of the reference and query structures.
			// and finally store the data in either vClaRNATuple or vSS;
			void readInputFiles(const fs::path& refPath, const fs::path& queryPath);
			
			
			
			// calcConfusionMatrix
			// gets two string for reference and query binary interaface interaction
			// and  calculate the confusion matrix's category.
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrix();
			
			// caclSimilarityScores gets a ConfusionMatrixTuple
			// and calculate all the score in the SimilarityScores class
			// and returns them as a map
			// in which the keys are the name of the scores and values are their value 
			[[nodiscard]]std::map<std::string, double> calcSimilarityScores(ConfusionMatrixTuple cmt);
			
			// wirte the results into a file
			void writeScores(ConfusionMatrixTuple cmt, fs::path outputPath);
			
			CompareStructures() = default;
			virtual ~CompareStructures() = default;
		
	};
}
