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
//#include "lib/SimilarityScores.h"
#include "lib/CompareStructures.h"

using ConfusionMatrixTuple_blast = std::tuple<int, int, int, int, int, int, int, int>; // (TP, TN, FP, FN, TPx, TNx, FPx, FNx);


struct BlastSSMap
{
	SSMap all;
	SSMap cans;
	SSMap nonCans;
};

namespace BLAST
{
	class CompareStructures : public SS::CompareStructures
	{
		public:
			std::string requestedInteractions;
			FindInteractions findInteractions;
			SimilarityScores simScores;
			BlastSSMap ref, query;
			
			
			// readStructures(blast)
			// by using the functions in the FindInteraction function read corresponding file
			// it gets paths of the reference and query structures.
			// and finally store the data in either vClaRNATuple or vSS;
			// it is needed to be called to initialized required variables in other functions
			void readFile(std::filesystem::path blastPath);
			
			// calcConfusionMatrix
			// gets the paths of the  blast file
			// by using readFile function will read structures and will compare the base-pairs
			[[nodiscard]] ConfusionMatrixTuple_blast calcConfusionMatrix();					
			[[nodiscard]] ConfusionMatrixTuple_blast calcConfusionMatrixCans();					
			[[nodiscard]] ConfusionMatrixTuple_blast calcConfusionMatrixNonCans();

			
			// caclSimilarityScores gets a ConfusionMatrixTuple
			// and calculate all the score in the SimilarityScores class
			// and returns them as a map
			// in which the keys are the name of the scores and values are their value 
			[[nodiscard]]std::map<std::string, double> calcSimilarityScores(ConfusionMatrixTuple_blast cmt, int max_n_positives);
			[[nodiscard]]std::vector<std::map<std::string, double>> calcSimilarityScores(std::vector<ConfusionMatrixTuple_blast> vcmt);
			
			
			std::vector<ConfusionMatrixTuple_blast> calcConfusionMatrix(fs::path path);
			
			// wirte the results into a file
			void writeScores(ConfusionMatrixTuple_blast cmt, fs::path outputPath);
			void writeScores(std::vector<ConfusionMatrixTuple_blast> vcmt, fs::path outputPath);
			
			CompareStructures();
			CompareStructures(std::string requestedIntercations);
			virtual ~CompareStructures() = default;
		private:
			std::vector<int> m_vMax_n_positives;
			int find_max_n_positives(const SSMap& ref, const SSMap& query);
		
	};
}


