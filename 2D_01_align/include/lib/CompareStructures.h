/* 
 * Author:  Seyed Naeim Moafinejad
 * Based on previous version by: Iswarya Pandara Nayaka PJ 
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
#include <optional>
#include "SimilarityScores.h"

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


namespace SS
{
	class CompareStructures
	{
		public:
			std::string requestedInteractions;
			std::vector<std::string> sequence;
			FindInteractions findInteractionRef, findInteractionQuery;
			SimilarityScores simScores;
			SSMap  vRef; // for storing extracted data from dot-bracket notation
			SSMap vQuery; // SSMap is a defined variable in FindInteraction.h
			
			// readStructures(
			// by using the functions in the FindInteraction class read corresponding file
			// it gets paths of the reference and query structures.
			// and finally store the data in either vClaRNATuple or vSS;
			void readStructures(fs::path refPath, fs::path queryPath);
			
			
			// calcConfusionMatrix
			// gets the paths of the  reference and query structures
			// by using readStructures function will read structures and will compare the base-pairs
			// for ss--> if a given base-pair was the same in both structures it would be TP
			//			 if there is no base-pair in the same place it would be TN
			//			 if there was a given base-pair in reference but not in query it would be FN
			//			 if there was not a given base-pair in reference but it was in query it would be FP							
			[[nodiscard]]ConfusionMatrixTuple calcConfusionMatrix(fs::path refPath, fs::path queryPath);
			[[nodiscard]]std::vector<ConfusionMatrixTuple> calcConfusionMatrix(fs::path refPath, fs::path queryPath, fs::path seqPath);
			[[nodiscard]]ConfusionMatrixTuple calcConfusionMatrixCans();
			[[nodiscard]]ConfusionMatrixTuple calcConfusionMatrixNonCans();
			
			// caclSimilarityScores gets a ConfusionMatrixTuple
			// and calculate all the score in the SimilarityScores class
			// and returns them as a map
			// in which the keys are the name of the scores and values are their value 
			std::map<std::string, double> calcSimilarityScores(ConfusionMatrixTuple cmt);
			
			// wirte the results into a file
			void writeScores(ConfusionMatrixTuple cmt, fs::path outputPath);
			void writeScores(std::vector<ConfusionMatrixTuple> vcmt, fs::path outputPath);
			
			CompareStructures();
			CompareStructures(bool withSeq, bool isWobble_canonical);
			virtual ~CompareStructures() = default;
		
	};
}
