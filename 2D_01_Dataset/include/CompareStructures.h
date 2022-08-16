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

using ScoreMatrix = std::vector<std::vector<double>>;
using ScoreMap = std::map<std::string, ScoreMatrix>;

namespace SS_all
{
	class CompareStructures : public SS::CompareStructures
	{
		
		public:
			std::string m_requestedScores;
			FindInteractions findInteractions;
			SimilarityScores simScores;
			ScoreMatrix scoreMatrix;
						
			//SSMap refMap, queryMap;
			//SSMatrix refMatrix, queryMatrix;
			
			
			bool get_is_matrix_prefered() const;
			void set_is_matrix_prefered(bool is_matrix_prefered);
			
			// readStructures(blast)
			// by using the functions in the FindInteraction function read corresponding file
			// it gets paths of the reference and query structures.
			// and finally store the data in vSS;
			// it is needed to be called to initialized required variables in other functions
			void readFile(const std::filesystem::path& allPath); // it think since we only have one input file or folder the reading file in the FindInteractions is enough
			
			// calcConfusionMatrix
			// gets SSMap or SSMatrix for reference and query structures
			// and  calculate the confusion matrix's category.
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrix(SSMap& ssMapRef, SSMap& ssMapQuery);
			[[nodiscard]] ConfusionMatrixTuple calcConfusionMatrix(const SSMatrix& ssMatrixRef, const SSMatrix& ssMatrixQuery);
			
			//calcSimilarityScores using the vector of different stuctures 
			// and calculate the requested similarity score.
			// and save the requested score in the matrix 
			// which shows requested score for structre i and j
			// the calculation is depend on the user 
			// if they prefer matrix calculation or vector calculation
			[[nodiscard]] ScoreMap calcSimilarityScoresVector();
			[[nodiscard]] ScoreMap calcSimilarityScoresMatrix();
			[[nodiscard]] ScoreMap calcSimilarityScores();
			
			
			//std::vector<ConfusionMatrixTuple_blast> calcConfusionMatrix(fs::path path);
			
			//// wirte the results into a file
			void writeScores(const std::filesystem::path& path);
			//void writeScores(std::vector<ConfusionMatrixTuple_blast> vcmt, fs::path outputPath);
			
			CompareStructures();
			CompareStructures(std::string requestedIntercations, bool is_matrix_prefered);
			virtual ~CompareStructures() = default;
	
		private:
			bool m_is_matrix_prefered;
			int m_max_n_positives = 0;
			void find_max_n_positives_vector();
			void find_max_n_positives_matrix();
		
	};
}

