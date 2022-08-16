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
#include <string>
#include <vector>


 
namespace Binary
{
	class FindInteractions
	{
		public:
			
			std::string binaryInterface;
			std::vector<std::string> m_vChains_binary;
			void readInputFile(const std::filesystem::path& path);
			int m_number_positives = 0;
			
			FindInteractions()=default;
			virtual ~FindInteractions()=default;
			
	};
}

