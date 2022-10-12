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
 
 /* import  required libraries */
 // include <format> // for time being it is only supported in the microsoft visual studio compiler
 // we can use format.dev. The foramt in standard library is a sublibrary of this external library

#include <cstddef>
#include <iostream>
#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

using SSMap = std::map<int, int>;
using SSContacts = std::vector<std::pair<int, int>>;
using SSMatrix = std::vector<std::vector<int>>;

//============================================================
// !!! Be careful with these two following types !!!         |
// each element of BasePairSymbols has BasePairSymbol type   | 
using BasePairSymbols = std::map<char, char>; //             |
using BasePairSymbol = std::pair<char, char>;//              |
//===============================================================================================================================================================================
										   //ChainID,NuclNumber,NuclName,         ChainID, NuclNumber, NuclName																    |
using ClaRNATuple = std::tuple</*first*/std::string, int, std::string, /*second*/std::string, int, std::string, /*type of interaction*/std::string, /*wieght*/ double>;//       |
										   // chainID   nuclNumber nucleName edge                   chainID,  nuclNumber nucleName  edge         cis or trans                   |
using FullClaRNATuple = std::tuple</*first*/std::string, int, std::string, std::string, /*second*/ std::string, int, std::string, std::string, std::string,/*weight*/ double>;//|
//===============================================================================================================================================================================
const std::string closeSymbols = ")]}>abcdefghijklmnopqrstuvwxyz";	  //|
const std::string openSymbols = "([{<ABCDEFGHIGKLMNOPQRSTUVWXYZ";	  //|
//=======================================================================	

// a structure for separating different symbols in in dot-bracket notation
struct OpenDotCloseIndices
{
	std::vector<int> vOpens, vDots, vCloses; 
};
 
namespace SS
{
	class FindInteractions
	{
		public:
			size_t sequence_length;
			size_t ssSize;
				
			// different  map for differnt types of interaction
			SSMap m_SSMap_can;
			SSMap m_SSMap_wobble;
			SSMap m_SSMap_noncan;
			SSMap m_SSMap;
			
			
			// if user choose to calculate the with contact matrix
			// the program use the following matrix;
			SSMatrix m_Matrix_can;
			SSMatrix m_Matrix_wobble;
			SSMatrix m_Matrix_noncan;
			SSMatrix m_Matrix_all;
			//SSMatrix ss_line; // it would be parsed into the matrix of different type
			
			// variables which are the same for both type of calculations
			std::string m_sequence;
			std::string m_seqWithSeparateChains;
			
			// setter and getter for private data
			int get_numberOfCanBasePairs() const;
			int get_numberOfNonCanBasePairs() const;
			int get_numberOfAllBasePairs() const;
			bool get_isWobble_canonical() const;
			void set_isWobble_canonical(bool isWobble_canonical) ;
			bool get_withSeq() const;
			void set_withSeq(bool withSeq);
			bool get_is_2D_on() const;\
			void set_is_2D_on(bool is_2D_on);
			
			bool isCanonical(std::pair<int, int> ssPair) const;
			bool isWobble(std::pair<int, int> ssPair) const;
			
			
			// read file and store the sequence and structures
			// all the structures have the same sequence
			void init_sequence(const std::filesystem::path& path);
			//initialized structures from single file
			void init_structure(const std::filesystem::path& path);			
			void readInputFiles(const std::filesystem::path& path);
			
			// it is not the override function in super class 
			// it gets the BasePairSymbols and vector of the secondary structures
			// separate the indices for the open and close brackets
			// at the end remains would be dots 
			[[nodiscard]] OpenDotCloseIndices separateOpensDotsCloses(const BasePairSymbols& basePairSymbols, const std::vector<std::string>& v_ss) const; 
			
			// ss2IndicesPair here is not the overriden here
			// it is the new one with different argument
			[[nodiscard]] SSContacts ss2IndicesPair(BasePairSymbols& basePairSymbols, std::vector<std::string> v_ss);
			
			[[nodiscard]] SSMap SSContacts2SSMap(const SSContacts& ssc);
			[[nodiscard]] SSMatrix SSContacts2SSMatrix(const SSContacts& ssc);
			
			// extractInteractions is different for the one in the SS::FindInteractions
			SSMap extractInteractions_vector(const std::vector<std::string>& v_ss);
			SSMatrix extractInteractions_matrix(const std::vector<std::string>& v_ss);
			
			
			// getting the path of the input and then read it and parse the informations
			void readInputFile(const std::filesystem::path& seqpath);
			
			// a fucntion for the separating canonical and noncanonical interactions
			void mapSS2seq(SSMap ssMap);
			
			
			
			// defualt costructor initialize public data member of the class
			// m_numberOfCanBasePairs = 0, m_numberOfNonCanBasePairs = 0, m_numberOfAllBasePairs = 0
			FindInteractions();
			
			
			// initializer constructor, initialize private data of the class by the user-defined value
			// and initialize public data member as in default c-tor
			FindInteractions(bool isWobble_canonical, bool m_withSeq);
		private:
			int m_numberOfCanBasePairs; 
			int m_numberOfNonCanBasePairs;
			int m_numberOfAllBasePairs;
			bool m_isWobble_canonical;
			bool m_withSeq;
			bool m_is_2D_on;
			
	};
}

void initSSMatrix(SSMatrix& ssMat, size_t length);

namespace stdExt{
	// mismatch only finds first mismatch and ignores others
	// so in the bellow you can find the extended version
	// which saves all the mismatches in a vector.
	template<typename InputIterator1, typename InputIterator2> 
	std::vector<std::pair<InputIterator1, InputIterator2>> mismatches(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2)
	{
		
		std::vector<std::pair<InputIterator1, InputIterator2>> vMismatches;
		while(first1 != last1 && first2 != last2)
		{
			if(*first1 != *first2)
			{
				vMismatches.push_back(std::make_pair(first1, first2));
			}
			++first1;
			++first2;
		}
		return vMismatches;
	}
	
	template<typename InputIterator1, typename InputIterator2, typename BinaryPredicate> 
	std::vector<std::pair<InputIterator1, InputIterator2>> mismatches(InputIterator1 first1, InputIterator1 last1, InputIterator2 first2, InputIterator2 last2, BinaryPredicate pred)
	{
		std::vector<std::pair<InputIterator1, InputIterator2>> vMismatches;
		while(first1 != last1 && first2 != last2)
		{
			if(pred(*first1,*first2))
			{
				vMismatches.push_back(std::make_pair(first1, first2));
			}
			++first1;
			++first2;
		}
		return vMismatches;
		
	}
	
	// find_if only finds the first match and igonres the others
	// this fuction is a extention of the find_if algorithm in std lib
	// and finds all the matches and inserts them in a vector
	template<typename InputIterator, typename Predicate>
	std::vector<InputIterator> find_all_if(InputIterator first, InputIterator last, Predicate pred)
	{
		std::vector<InputIterator> dest;
		while(first != last)
		{
			if(pred(*first))
			{
				dest.push_back(first);
			}
			++first;
		}
		
		return dest;
	}
	
	// a template helper function and functions for printing tuples
	template<typename Tuple, size_t... Indices>
	void tuplePrint(const Tuple& t, std::index_sequence<Indices...>)
	{
		((std::cout << std::get<Indices>(t) << " "), ...);
	}
	
	template<typename... Args>
	void tuplePrint(const std::tuple<Args...>& t)
	{
		tuplePrint(t, std::index_sequence_for<Args...>());
		std::cout << std::endl;
	}
}

