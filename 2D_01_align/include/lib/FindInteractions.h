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
			
			int ssSize;
			
			std::string seqWithSeparateChains;
			
			// diferent  map for differnt types of interaction
			SSMap ssMap_can;
			SSMap ssMap_noncan;
			
			// setter and getter for private data
			int get_numberOfCanBasePairs() const;
			int get_numberOfNonCanBasePairs() const;
			int get_numberOfAllBasePairs() const;
			 
			bool get_isWobble_canonical() const;
			void set_isWobble_canonical(bool isWobble_canonical) ;
			bool get_withSeq() const;
			void set_withSeq(bool withSeq);
			
			// isCanonical is helper function
			// to find canonical interaction 
			// with mapping ssPair on sequence
			bool isCanonical(std::vector<std::string> vSeqOfChains, std::pair<int, int> ssPair) const;
			bool isCanonical(std::string Seq, std::pair<int, int> ssPair) const;
			
			// separateOpensDotsCloses is a helper function 
			// that gets BasePairSymbols and secondary structure string
			// separate different Symbols indices and returns an structure which contains all of them
			// it throw exception it would be better to use with try and catch
			[[nodiscard]]OpenDotCloseIndices separateOpensDotsCloses(BasePairSymbols basePairSymbols, std::string ss) const;
			
			// ss2IndicesPair is a helper function
			// get OpenDotCloseIndices 
			// and use the information saved in this structure 
			// to make pair of indices that are involved in a base pairs
			// if a nucleotide does not have contribution in base pairtin
			// it will take -1 as its couple in base pair
			// and it returns SSMap
			// it throws exception it would be better to use it with try and catch.
			[[nodiscard]]SSMap ss2IndicesPair(BasePairSymbols basePairSymbols, std::string ss); 
			
			// extract interaction gets a string of secondary structure
			// and parse its informations and save them into a SSMap
			[[nodiscard]]SSMap extractInteractions(std::string ss);
			
			// following function overloads extract interaction
			// to work for multiple-line secondary structure;
			// and returns a vector of SSMap,
			// information of each line will be save in different element;
			// I'm not sure to use it or not ???
			[[nodiscard]]std::vector<SSMap> extractInteractions(std::vector<std::string> ss);
			
			// readInputFile get a path and read secondary structure and returns SSMap (for time being)
			[[nodiscard]]SSMap readInputFile(std::filesystem::path path);
			[[nodiscard]]std::vector<std::string> readSeqFile(std::filesystem::path path);
			
			// mapSS2Seq map structure into sequence
			// and separate the canonical and non-canonical interaction 
			// it gets sequence and ssMap
			void mapSS2seq(std::vector<std::string> vSeqOfChains, SSMap ssMap);
			void mapSS2seq(std::string seq, SSMap ssMap);
			
			// defualt costructor initialize public data member of the class
			// m_numberOfCanBasePairs = 0, m_numberOfNonCanBasePairs = 0, m_numberOfAllBasePairs = 0
			FindInteractions();
			
			// initializer constructor, initialize private data of the class by the user-defined value
			// and initialize public data member as in default c-tor
			FindInteractions(bool isWobble_canonical, bool withSeq);
			
			
			
			
		private:		
			// numberOfCanBasePairs would be extracted from both type of inputs (ss and clarna) with/without sequence
			int m_numberOfCanBasePairs; 
			// if sequence for the structure is available it could be extracted from both type of inputs
			// otherwise it only can be extracted from clarna outputs as inputs
			int m_numberOfNonCanBasePairs;
			// numberOfAllInteracions extracted from all kind of inputs with/without sequence
			int m_numberOfAllBasePairs;
		
			// variables that indiacte if the 
			bool m_isWobble_canonical = false;
			bool m_withSeq = false;
			
	};
}

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

