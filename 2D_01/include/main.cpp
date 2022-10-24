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


#include <algorithm>
#include "CompareStructures.h"
#include <getopt.h>
#include <chrono>

using namespace std;
using namespace std::chrono;


void showHelp(){
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore_2D_01\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-r or --reference: reference file\n";
	cout << "-q or --query: target file (it can be bracket dot notation of SS or ClaRNA.out depending on option )\n";
	//cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	cout << "-d or --dotBracket: for comparison *.ss files with each other\n";
	cout << "\t (the type of the interaction must be chosen)\n";
	cout << "\t (c = canonical or e = canonical + wobble (extended), w = wobble, n = non-canonical, a = all)\n";
	cout << "\t (C = canonical or E = canonical + wobble (extended), W = wobble, N = non-canonical, A = all)\n";
	cout << "\t you cannot ask for canonical (c) and extended canonical (canonical +  wobble) at the same time\n";
	cout << "\t you can choose more than one interaction for example -d cna";
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-v or --1D: comparison with 1D array algorithm\n";
	cout << "-m or --2D: comparison with 2D array algorithm";
	cout << "-s or --sequence: sequence file\n";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	

	cout << "Usage:\n\t./2D_01 -r <referenceFile> -s [sequenceFile] -q <queryFile> -d <requested_interactions> --1D (or --2D) -o [outputName]";
	cout << "\nexample:\n\t./2D_01  -r samples/dotBracketRef.SS -q samples/dotBracketQuery.SS -s samples/SeqForDotBracket.seq -d enaw --1D -o results/sampleTest.csv";
	cout << endl;
	exit(EXIT_FAILURE);
}

void printRequestedInteractions(string requestedInteractions)
{
	auto itC = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'C' || ch == 'c');});
	auto itE = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' || ch == 'e');});
	
	if(itC != cend(requestedInteractions) && itE != cend(requestedInteractions))
	{
		cout << " you cannot ask for canonical (c) and extended canonical (canonical +  wobble) at the same time\n"; 
		exit(EXIT_FAILURE);
	}
	
	if(itC != cend(requestedInteractions))
	{
		cout << "\tCanonical" << endl;
	}
	
	if(itE != cend(requestedInteractions))
	{
		cout << "\tCanonical + wobble (extended canonical)" << endl;
	}
		
	auto itW = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'W' || ch == 'w');});
	if(itW != cend(requestedInteractions))
	{
		cout << "\tWobble" << endl;
	}
	
	auto itN = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'N' || ch == 'n');});
	if(itN != cend(requestedInteractions))
	{
		cout << "\tNon-canonical" << endl;
	}
	
	auto itA = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'A' || ch == 'a');});
	if(itA != cend(requestedInteractions))
	{
		cout << "\tAll" << endl;
	}
}


int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
		exit(EXIT_FAILURE);
	}
	
	
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore_2D_01\n";
	cout << "===============================\n";
	
    string refPath, queryPath, seqFile, outFilename; 
    bool dotBracketFlag, ClaRNAFlag, outFilenameFlag, refFlag, queryFlag, sequenceFlag, vectorFlag, matrixFlag;
    dotBracketFlag = ClaRNAFlag = outFilenameFlag = refFlag = queryFlag = sequenceFlag = vectorFlag = matrixFlag = false;

    const struct option long_opts[] = {
    	{"query", required_argument, nullptr,'q'},
    	{"reference", required_argument, nullptr, 'r' },
    	{"sequence", required_argument, nullptr, 's'},
    	{"dotBracket", no_argument,  nullptr, 'd'},
    	{"1D", no_argument, nullptr, 'v'},
    	{"2D", no_argument, nullptr, 'm'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	string requestedInteractions = "acn";
	int optind = 0;

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"q:r:s:d:o:vmh", long_opts, &optind);
        if (-1 == opt)
        {
            break;
		}
        
        switch (opt)
        {
			case 'q':
				queryFlag = true;
				queryPath = optarg; 
				cout << "Query path: " << optarg << endl;
				break;
			case 's':
				sequenceFlag = true;
				seqFile=optarg; 
				cout << "Input sequence file: " << seqFile << endl;
				break;
			case 'r':
				refFlag = true;
				refPath = optarg;
				cout << "Reference path: " << optarg << endl;
				break;
				
			case 'd':
				dotBracketFlag = true;
				requestedInteractions = optarg;
				cout << "Comparison of dot bracket notation will be done\n";
				cout << "Your requested interactions: \n";
				printRequestedInteractions(requestedInteractions);
				break;
			
			case 'v':
				vectorFlag = true;
				cout << "1D array algorithm is chosen" << endl;
				if(matrixFlag)
				{
					cout << "Error: 2D algorithm is already chosen, you cannot choose both algorithm\n";
					exit(EXIT_FAILURE);
				}
				break;
			
			case 'm':
				matrixFlag = true;
				cout << "2D array algorithm is chosen" << endl;
				if(vectorFlag)
				{
					cout << "Error:  1D algorithm is already chosen, you cannot choose both algorithm\n";
					exit(EXIT_FAILURE);
				}
				break;
			
			case 'o':
				outFilenameFlag = true;
				outFilename = optarg;
				break;
			case 'h': 
			case '?': 
			default:
				showHelp();
				break;
		}
    }
    
    // cheking input files
    if(!refFlag)
    { 
		cout << "Reference file is required\n";
		showHelp();
	}
	
	if(!queryFlag)
	{
		cout << "Query file is required\n";
		showHelp();
	}
	
	
	if(!outFilenameFlag)
	{
		outFilename = "FinalResuts.csv";
	}
	
	auto start = chrono::high_resolution_clock::now();
	if(dotBracketFlag)
	{
		SS::CompareStructures cs;
		cs.m_requestedInteractions = requestedInteractions;
		auto itE = find_if(requestedInteractions.begin(), requestedInteractions.end(), [](char ch){return(ch == 'E' || ch == 'e');});
		if(itE != requestedInteractions.end())
		{
			cs.findInteractionRef.set_isWobble_canonical(true);
			cs.findInteractionQuery.set_isWobble_canonical(true);
		}
		
		cs.set_is_2D_on(matrixFlag);
		cs.readStructures(refPath, queryPath);
		//cs.findInteractionRef.set_withSeq(sequenceFlag);
		//cs.findInteractionQuery.set_withSeq(sequenceFlag);
		std::vector<ConfusionMatrixTuple> vcmt;
		
		if(!sequenceFlag)
		{
			if(requestedInteractions != "a")
			{
				cout << "Without sequence file only 'a or A' request is available" << endl;
				cs.m_requestedInteractions = "a";
				//exit(EXIT_FAILURE);
			}
			
			if(vectorFlag)
			{
				vcmt = cs.calcConfusionMatrixVector();
			}
			
			if(matrixFlag)
			{
				vcmt = cs.calcConfusionMatrixMatrix();
			}
		}
		else
		{
			cs.readsequence(seqFile);
			if(vectorFlag)
			{
				vcmt = cs.calcConfusionMatrixVector();
			}
			
			if(matrixFlag)
			{
				vcmt = cs.calcConfusionMatrixMatrix();		
			}
		}
		cs.writeScores(vcmt, outFilename);
	}
	else
	{
		cout << "Please add opton -d in your command line" << endl;
		cout << "This optoin is obligatory for integrity of usage whole software" << endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	
	auto end = high_resolution_clock::now();
	auto duration_s = duration_cast<seconds>(end - start);
	auto duration_ms = duration_cast<milliseconds>(end - start);
	
	if(duration_s.count() == 0) cout << "running time: " << duration_ms.count() << " ms\n";
	else cout << "running time: " << duration_s.count() << " s\n";
	
	cout << "Done ;)" << endl;

	
	
	return 0;
}
