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
#include <filesystem>
#include <iostream>
#include <iterator>
#include "ReadPDB.h"
#include <string>
#include <string_view>
#include <utility>
#include <vector>

//===============================================================================================================================================================================
										   //ChainID,NuclNumber,NuclName,         ChainID, NuclNumber, NuclName																    |
using ClaRNATuple = std::tuple</*first*/std::string, int, std::string, /*second*/std::string, int, std::string, /*type of interaction*/std::string, /*wieght*/ double>;//       |
										   // chainID   nuclNumber nucleName edge                   chainID,  nuclNumber nucleName  edge         cis or trans                   |
using FullClaRNATuple = std::tuple</*first*/std::string, int, std::string, std::string, /*second*/ std::string, int, std::string, std::string, std::string,/*weight*/ double>;//|
//===============================================================================================================================================================================


using std::cout;
using std::cerr;
using std::endl;

const int LEAST_POSSIBLE_INTERACTIONS_CANS     = 1;
const int LEAST_POSSIBLE_INTERACTIONS_WOBBLES  = 1;
const int LEAST_POSSIBLE_INTERACTIONS_NONCANS  = 3;
const int LEAST_POSSIBLE_INTERACTIONS_STACKING = 2; // also 3 is impossible for stackin
												    // and with other type interactions only 5 is acceptable, it would be checked in the functions.

const int ONE_EDGE    = 1; // only for ww_cis (canonical and wobble)
const int TWO_EDGES   = 2; // onlly for stackings
const int THREE_EDGES = 3; // for canonical, noncanonical, wobbles, all base pairs, and any combinations of them
const int FIVE_EDGES  = 5; // for all interactions. 


const int WW_CIS  	=  1;
const int WW_TRAN 	=  2;
const int WH_CIS  	=  3;
const int WH_TRAN 	=  4;
const int WS_CIS  	=  5;
const int WS_TRAN 	=  6;

const int HW_CIS  	=  7;
const int HW_TRAN 	=  8;
const int HH_CIS  	=  9;
const int HH_TRAN 	= 10;
const int HS_CIS  	= 11;
const int HS_TRAN   = 12;

const int SW_CIS    = 13;
const int SW_TRAN   = 14;
const int SH_CIS    = 15;
const int SH_TRAN   = 16;
const int SS_CIS    = 17;
const int SS_TRAN   = 18;

const int STACK_3_3 = 19;
const int STACK_3_5 = 20;
const int STACK_5_3 = 21;
const int STACK_5_5 = 22;
					
//	   						   Wc	Wt	 Hc   Ht   Sc   St    
using Edge_tuple = std::tuple<int, int, int, int, int, int>;
//                             >    <
using Face_tuple = std::tuple<int, int>;
//										W			 H			  S		  	>            < 
using ClaRNAMatrix_elem = std::tuple<Edge_tuple, Edge_tuple, Edge_tuple, Face_tuple, Face_tuple>;
using ClaRNAMatrix = std::vector<std::vector<ClaRNAMatrix_elem>>;


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
		((cout << std::get<Indices>(t) << " "), ...);
	}
	
	template<typename... Args>
	void tuplePrint(const std::tuple<Args...>& t)
	{
		tuplePrint(t, std::index_sequence_for<Args...>());
		cout << endl;
	}
	
	template<typename... Args>
	void tuplePrint2(const std::tuple<Args...>& t)
	{
		tuplePrint(t, std::index_sequence_for<Args...>());
		cout << "  ";
	}
}

namespace CLARNA
{
	class FindInteractions
	{
		public:	 
			std::string seqWithSeparateChains;
			std::string m_sequence;
			
			// vectors of different interactions which would be initilized by separate functions
			std::vector<ClaRNATuple> vtCans;
			std::vector<ClaRNATuple> vtWobbles;
			std::vector<ClaRNATuple> vtStacks;
			std::vector<ClaRNATuple> vtNonCans;
			std::vector<ClaRNATuple> vtBasePairs;
			std::vector<FullClaRNATuple> vftCans;
			std::vector<FullClaRNATuple> vftWobbles;
			std::vector<FullClaRNATuple> vftStacks;
			std::vector<FullClaRNATuple> vftNonCans;
			std::vector<FullClaRNATuple> vftBasePairs;
			std::vector<FullClaRNATuple> vftAll;
			
			// matrices of different interactions which would be initialized by separate functinos
			ClaRNAMatrix m_Matrix_can;
			ClaRNAMatrix m_Matrix_wobble;
			ClaRNAMatrix m_Matrix_noncan;
			ClaRNAMatrix m_Matrix_baseParis;
			ClaRNAMatrix m_Matrix_stacks;
			ClaRNAMatrix m_Matrix_all;
			
			std::vector<Atom> m_vNucleotides;
			
			// setter and getter for private data
			int get_numberOfCanBasePairs() const;
			int get_numberOfNonCanBasePairs() const;
			int get_numberOfAllBasePairs() const;
			int get_numberOfStacks() const;
			int get_numberOfAllInteractions() const;
			 
			bool get_isWobble_canonical() const;
			void set_isWobble_canonical(bool isWobble_canonical);
			
			int get_number_involved_faces_edges() const;
			void set_number_involved_faces_edges(int number_involved_faces_edges);
			
			bool isStacking(const ClaRNATuple& ct) const;
			bool isStacking(std::string nuclEdge) const;
			bool isCanonical(const ClaRNATuple& ct)const;
			bool isCanonical(const FullClaRNATuple& fct) const;
			bool isWobble(const ClaRNATuple& ct) const;
			bool isWobble(const FullClaRNATuple& fct) const;
			bool isBasePhosphate(const ClaRNATuple& fct) const;
			bool isBaseRibose(const ClaRNATuple& fct) const;
			bool isBaseBackbone(const ClaRNATuple& ct) const;
	
			
			// extractInteraction gets a line of clarna out
			// and parse its informations and save them into a ClaRNATuple
			[[nodiscard]]ClaRNATuple extractInteractions(std::string clarnaLine);
			
			// readInputFile gets a path and read interactions detected by clarna and returns a vector of ClaRNATuple
			[[nodiscard]]std::vector<ClaRNATuple> readInputFile(std::filesystem::path path);
			
			// readPDBFile gets a path and extract pdb information and return a fully non-intractde FullClaRNATuple
			// later we can use the FullClaRNATuple to compare the FullClaRNATuple extracted from clarna
			// if the were not the same the program can throw an error
			[[nodiscard]]std::vector<FullClaRNATuple> readPDBFile(std::filesystem::path path);
			
			// makeFullyNonInteractedVector read pdb file and make fully noninteracted vector
			// for all the nucleotide in different modes
			// when user asks for the presentes of all 5 edges
			[[nodiscard]]std::vector<FullClaRNATuple> makeFullyNonInteractedVector5();
			// when user ask for only 1 edge (possible for canonical interactions)
			[[nodiscard]]std::vector<FullClaRNATuple> makeFullyNonInteractedVector1();
			// when user ask for 3 edges
			[[nodiscard]]std::vector<FullClaRNATuple> makeFullyNonInteractedVector3();
			//when user ask for 2 edges (possible for stacking)
			[[nodiscard]]std::vector<FullClaRNATuple> makeFullyNonInteractedVector2();
			
						
			// a function for separating canonical base pairs from others
			void sepCanBasePairs(const std::vector<ClaRNATuple>& vtInteractions, const std::vector<FullClaRNATuple>& vfct_formPDB);
			//// the difference with previous function is, all other interactions are considered as TN.
			//void sepCanBasePairs_TN_All(const std::vector<ClaRNATuple>& vtInteractions, const std::vector<FullClaRNATuple>& vfct_formPDB);
			//// the difference with previous function is, noncanonical interactions are considered as TN.
			//void sepCanBasePairs_TN_nonCans(const std::vector<ClaRNATuple>& vtInteractions, const std::vector<FullClaRNATuple>& vfct_formPDB);
			//// the difference with previous function is, stacking interactions are considered as TN.
			//void sepCanBasePairs_TN_stacking(const std::vector<ClaRNATuple>& vtInteractions, const std::vector<FullClaRNATuple>& vfct_formPDB);
			
			
			
			// a function for separating non-canonical base pairs from others
			void sepNonCanBasePairs(const std::vector<ClaRNATuple>& vtInteractions, const std::vector<FullClaRNATuple>& vfct_formPDBs);
			
			// a function for separating stacking from other type of interatctions
			void sepStacking(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB);
			
			// a function for separating all type of the base-pairs from stacking
			void sepBasePairs(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB);
			
			// a function for separating wobble interactions.
			void sepWobble(const std::vector<ClaRNATuple>& vct, const std::vector<FullClaRNATuple>& vfct_formPDB);
			
			// a function for comparing two clarna tuple
			// it would be use for comparintg the FullClaRNATuples from clarna and pdb
			void checkClaRNAwithPDB(const std::vector<FullClaRNATuple>& vfct_clarna, const std::vector<FullClaRNATuple>& vfct_pdb);
			
				
			// this function gets the vector of ClaRNATuple 
			// and parses the edge from type of interaction for each nucleotide involved in the base-base interactions;
			[[nodiscard]]std::vector<FullClaRNATuple> makeFullVectorOfClaRNATuple(const std::vector<ClaRNATuple>& vct);
				
			// makeFullClaRNATuple is A helper function to get a ClaRNATuple 
			// and returns FullClaRNATuple 
			FullClaRNATuple  makeFullClaRNATuple(const ClaRNATuple& ct);
			
			// makeNonInteractedVector gets a vector of the sequenec of each chain
			// and makes a vector of FullClaRNATuple without any type of base inteactions
			// and returns that vector
			std::vector<FullClaRNATuple> makeNonInteractedVector(const std::vector<std::string>& vSeq)const;
				
			// stackingReciprocal gets a string
			// checks the type of the stacking for the nucleotide
			// return the opposite type of stacking (example ">" --> "<")
			std::string stackingReciprocal(std::string nuclEdge) const;
			
			// makeReciprocalInteractions gets a vector of FullClaRNATuple 
			// and make reciprocal inteaction out of any base interaction
			// and return new FullClaRNATuple
			[[nodiscard]]std::vector<FullClaRNATuple> makeReciprocalInteractions(const std::vector<FullClaRNATuple>& vfct) const;
			
			//make1Dvector
			
			[[nodiscard]]std::vector<FullClaRNATuple> make1Dvector(const std::vector<FullClaRNATuple>& vfct) const;
			
			// checkFirstNuclEdges gets a vector of FullClaRNATuple (or path to make it) and vector of sequence (or path to make it)
			// and checks the existing edge of each nucleotide which appeared as first nucleotide in the interactions.
			// returns a vector of BooleanEdgeTuple
			// since we have the function to make reciprocal interactions
			// only checking the first nucleotide is ok
			[[nodiscard]]std::vector<FullClaRNATuple> addNonInteractedEdges(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_clarna) const;
			[[nodiscard]]std::vector<FullClaRNATuple> addNonInteractedEdgesForCans(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_clarna) const;
			[[nodiscard]]std::vector<FullClaRNATuple> addNonInteractedEdgesForNonCans(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_clarna) const;
			[[nodiscard]]std::vector<FullClaRNATuple> addNonInteractedEdgesForStacks(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_clarna) const;
			[[nodiscard]]std::vector<FullClaRNATuple> addNonInteractedEdgesForBasePairs(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_clarna) const;
			[[nodiscard]]std::vector<FullClaRNATuple> addNonInteractedEdgesForWobble(const std::vector<FullClaRNATuple>& vfct_fromPDB, const std::vector<FullClaRNATuple>& vfct_clarna) const;
			
			void addNonInteractedEdges(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath);
			void addNonInteractedEdgesForCans(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath);
			void addNonInteractedEdgesForNonCans(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath);		
			void addNonInteractedEdgesForStacks(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath);
			void addNonInteractedEdgesForBasePairs(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath);
			void addNonInteractedEdgesForWobble(const std::filesystem::path& pdbPath, const std::filesystem::path& clarnaPath);
			
			// finalBasePairCheck gets a vector of FullClaRNATuple (with reciprocal interactions) 
			// and check the base pairs, if it found any base pair without reciprocal interaction (or vice versa)
			// it would be replace the base pair with non interacted nucleotide.
			// and finally it returns a vector of FullClaRNATuple
			[[nodiscard]]std::vector<FullClaRNATuple> finalBasePairCheck(const std::vector<FullClaRNATuple>& vfct);
			[[nodiscard]]std::vector<FullClaRNATuple> finalStackingCheck(const std::vector<FullClaRNATuple>& vfct);
			
			// overloaded finalBasePairCheck gets FullClaRNATuple additoiinally to check each base separately
			[[nodiscard]]FullClaRNATuple finalBasePairCheck(const std::vector<FullClaRNATuple>& vfct, const FullClaRNATuple& fct);
			[[nodiscard]]FullClaRNATuple finalStackingCheck(const std:: vector<FullClaRNATuple>& vfct, const FullClaRNATuple& fct);
			
			
			[[nodiscard]]ClaRNAMatrix v_ClaRNATuple2ClaRNAMatrix(const std::vector<ClaRNATuple>& vct) const;
			
			// defualt costructor initialize public data member of the class
			// m_numberOfCanBasePairs = 0, m_numberOfNonCanBasePairs = 0, m_numberOfAllBasePairs = 0, m_numberOfStacks = 0, m_numberOfAllInteractions = 0;
			FindInteractions();
			
			// initializer constructor, initialize private data of the class by the user-defined value
			// and initialize public data member as in default c-tor
			FindInteractions(bool isWobble_canonical);
			
			~FindInteractions() = default;
			
			
			
			
		private:
			int m_numberOfCanBasePairs; 
			int m_numberOfNonCanBasePairs;
			int m_numberOfAllBasePairs;
			// numberOfStacking extracted from only clarna ouptus
			int m_numberOfStacks;
			int m_numberOfAllInteractions;
			bool m_isWobble_canonical = false;
			int m_number_involved_faces_edges;
	};
	
}
