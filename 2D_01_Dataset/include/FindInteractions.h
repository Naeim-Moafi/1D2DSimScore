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

#include "lib/FindInteractions.h"
using SSContacts = std::vector<std::pair<int, int>>;
using SSMatrix = std::vector<std::vector<int>>;

namespace SS_all
{
	class FindInteractions : public SS::FindInteractions
	{
		public:
			size_t sequence_length;
				
			// different  map for differnt types of interaction
			std::vector<SSMap> m_vSSMaps_can;
			std::vector<SSMap> m_vSSMaps_noncan;
			std::vector<SSMap> m_vSSMaps;
						
			// if user choose to calculate the with contact matrix
			// the program use the following matrix;
			std::vector<SSMatrix> m_vMatrix_can;
			std::vector<SSMatrix> m_vMatrix_noncan;
			std::vector<SSMatrix> m_vMatrix_all;
			//SSMatrix ss_line; // it would be parsed into the matrix of different type
			
			// variables which are the same for both type of calculations
			std::string m_sequence;
			std::vector<std::string> m_vStructures_name;
			
			// setter and getter for private data
			int get_numberOfCanBasePairs() const;
			int get_numberOfNonCanBasePairs() const;
			int get_numberOfAllBasePairs() const;
			bool get_isWobble_canonical() const;
			void set_isWobble_canonical(bool isWobble_canonical);
			bool get_seqFlag() const{ return m_seqFlag;}
			std::string get_extension() const;
			void set_extension(const std::string& extension);
						
			// read file and store the sequence and structures
			// all the structures have the same sequence
			void init_sequence(const std::filesystem::path& path);
			//initialized structures from single file
			void init_structures_single(const std::filesystem::path& path);
			// initialized structures from multiple files
			void init_structures_multiple(const std::filesystem::path& path);
			
			void init_structures(const std::filesystem::path& path);
			
			// the path in this function is a file with several different structures
			void readSingleFile(const std::filesystem::path& path);
			
			// the path is a folder that contains several files for different structures
			void readMultipleFiles(const std::filesystem::path& path);
			
			// defualt costructor initialize public data member of the class
			// m_numberOfCanBasePairs = 0, m_numberOfNonCanBasePairs = 0, m_numberOfAllBasePairs = 0
			FindInteractions();
			
			// initializer constructor, initialize private data of the class by the user-defined value
			// and initialize public data member as in default c-tor
			FindInteractions(bool isWobble_canonical);
		private:
			int m_numberOfCanBasePairs;
			int m_numberOfNonCanBasePairs;
			int m_numberOfAllBasePairs;
			bool m_isWobble_canonical;
			bool m_seqFlag = false;
			bool m_extFlag = false;
			std::string m_extension;
			
	};
}


//#endif //__FINDINTERACTIONS_BLAST__
