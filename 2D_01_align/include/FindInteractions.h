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

 
 #ifndef __FINDINTERACTIONS_BLAST__
 #define __FINDINTERACTIONS_BLAST__
 
 /* import  required libraries */
 // include <format> // for time being it is only supported in the microsoft visual studio compiler
 // we can use format.dev. The foramt in standard library is a sublibrary of this external library

#include "lib/FindInteractions.h"

struct Range {
	int start;
	int term;
};

namespace BLAST
{
	class FindInteractions : public SS::FindInteractions
	{
		public:
						
			// different  map for differnt types of interaction
			std::vector<SSMap> m_vSSMaps_can;
			std::vector<SSMap> m_vSSMaps_noncan;
			std::vector<SSMap> m_vSSMaps;
			std::vector<std::string> m_vSequences;
			std::vector<std::string> m_vStructures;
			std::vector<std::string> m_vCansStructures;
			std::vector<std::string> m_vNonCansStructures;
			std::vector<Range> m_vRanges;
			
			// setter and getter for private data
			int get_numberOfCanBasePairs() const;
			int get_numberOfNonCanBasePairs() const;
			int get_numberOfAllBasePairs() const;
			bool get_isWobble_canonical() const;
			void set_isWobble_canonical(bool isWobble_canonical) ;
			
			void init_ranges(std::ifstream& inpFile);
			void init_structures(std::ifstream& inpFile);
			void init_sequences(std::ifstream& inpFile);
			int split_string(std::string_view string2parse, const char& delim, std::vector<std::string>& output_list) const;
			void readInputFile(std::filesystem::path path);
			[[nodiscard]]SSMap resetIndices(Range range, SSMap ssMap) const;
			[[nodiscard]]std::vector<SSMap> resetIndices(std::vector<Range> vRanges, std::vector<SSMap> vSSMap) const;
			std::string structure2Can(const std::string& sequence, const std::string& structure, SSMap ssMap) const;
			void structure2Can();
			std::string structure2NonCan(const std::string& sequence,const std::string& structure, SSMap ssMap) const;
			void structure2NonCan();
			void showInformation() const;
			void showInformationInRange(Range range, SSMap ssMap, int structure_number, const std::string& interaction_tye="A")const;
			
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
			int m_isWobble_canonical;
			
	};
}


#endif //__FINDINTERACTIONS_BLAST__
