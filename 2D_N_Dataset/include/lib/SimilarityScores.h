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


 
#ifndef _SIMILARITY_SCORES_H_
#define _SIMILARITY_SCORES_H_

// import required header and libraries
#include <iostream>



class SimilarityScores
{
	public:
		int TP, TN, FP, FN; // confusion matrix categories
		
		
		// setCategories
		// a function for setting TP, TN, FP, FN when an object is defined by default ctor
		void setCategories(int TP, int TN, int FP, int FN);
		
		/*** similarity scores ***/
		
		
		// mcc has range in the interval [-1, +1]
		// -1 is prefect misclassification 
		// +1 is prefect calssification
		// MCC produce high score if prediction obtained good results
		// in all four categories of confusion matrix (TP, TN, FP, FN)
		[[nodiscard]] double calcMCC() const;
		
		// Precistion or Positive Prediction Value(PPV) 
		// is a proportion of positive results in statistics
		// A high result can be interpreted as indicating the accuracy of such a statistic
		[[nodiscard]] double calcPrecision() const;
				
		// Recall or sensivity 
		// refers to the proportion of those that received a positive prediction on the test out of those that actually are positive
		[[nodiscard]] double calcRecall() const;
		
		// Fscore
		// it is the measure of test's accuracy
		// it would be calculated from precision and recall
		[[nodiscard]] double calcFscore() const;
		
		// FMIndex or Fowlkes–Mallows index (mean in previous version)
		// is an external evaluation method that is used to determine the similarity between two structures (here).
		// it also would be calculated from precision and recall 
		[[nodiscard]] double calcFMIndex() const;
		
		// JIndex
		[[nodiscard]] double calcJIndex() const;
		
		// Specificity
		[[nodiscard]] double calcSpecificty() const;
		
		// Balance Accuracy
		[[nodiscard]] double calcBA() const;
		
		// False omision rate
		[[nodiscard]] double calcFOR() const;
		
		// Prevalence threshold
		[[nodiscard]] double calcPT() const;
		
		// Critical success index
		[[nodiscard]] double calcCSI() const;
		
		// Markedness
		[[nodiscard]] double calcMK() const;
		
		SimilarityScores(){};
		SimilarityScores(int TP, int TN, int FP, int FN);
		//SimilarityScores(int&& TP, int&& TN, int&& FP, int&& FN);
		virtual ~SimilarityScores() = default;
	
	
};


#endif //_SIMILARITY_SCORES_H_
