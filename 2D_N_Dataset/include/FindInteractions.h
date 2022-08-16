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
#include "lib/FindInteractions.h"

namespace CLARNA_all
{
	class FindInteractions : public CLARNA::FindInteractions
	{
		public:	 
			std::string seqWithSeparateChains;
			std::string m_sequence;
			std::vector<std::string> m_vStructures_name;
			
			// vectors of different interactions which would be initilized by separate functions
			std::vector<std::vector<FullClaRNATuple>> v_vftCans;
			std::vector<std::vector<FullClaRNATuple>> v_vftWobbles;
			std::vector<std::vector<FullClaRNATuple>> v_vftStacks;
			std::vector<std::vector<FullClaRNATuple>> v_vftNonCans;
			std::vector<std::vector<FullClaRNATuple>> v_vftBasePairs;
			std::vector<std::vector<FullClaRNATuple>> v_vftAll;
			
			// matrices of different interactions which would be initialized by separate functinos
			//ClaRNAMatrix m_Matrix_can;
			//ClaRNAMatrix m_Matrix_wobble;
			//ClaRNAMatrix m_Matrix_noncan;
			//ClaRNAMatrix m_Matrix_baseParis;
			//ClaRNAMatrix m_Matrix_stacks;
			//ClaRNAMatrix m_Matrix_all;
			
			
			bool get_isWobble_canonical() const;
			void set_isWobble_canonical(bool isWobble_canonical);
			int get_number_involved_faces_edges() const;
			void set_number_involved_faces_edges(int number_involved_faces_edges);
			std::string get_extension() const;
			void set_extension(const std::string& extension);
			
			void readInputFiles(const std::filesystem::path& inputPath);

			
			// defualt costructor initialize public data member of the class
			// m_numberOfCanBasePairs = 0, m_numberOfNonCanBasePairs = 0, m_numberOfAllBasePairs = 0, m_numberOfStacks = 0, m_numberOfAllInteractions = 0;
			FindInteractions();
			
			// initializer constructor, initialize private data of the class by the user-defined value
			// and initialize public data member as in default c-tor
			FindInteractions(bool isWobble_canonical);
			
			~FindInteractions() = default;
			
			
			
			
		private:
			bool m_isWobble_canonical = false;
			int m_number_involved_faces_edges;
			std::string m_extension = ".out";
			bool m_pdbFlag;
			bool m_extFlag = false;
	};
	
}
