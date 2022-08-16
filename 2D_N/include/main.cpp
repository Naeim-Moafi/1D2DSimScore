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

using namespace std;


void showHelp(){
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore_2D_N\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-r or --reference: reference file (it can be bracket dot notation of SS or ClaRNA.out depending on option )\n";
	cout << "-q or --query: target file (it can be bracket dot notation of SS or ClaRNA.out depending on option )\n";
	cout << "-p or --pdb: sequence file (required for clarna results)\n";
	cout << "-c or --clarna_out:  for comparison *.out files with each other\n";
	cout << "\t (the type of the interaction must be chosen (there is no defualt type))\n";
	cout << "\t (c = canonical or e = canonical + wobble (extended), w = wobble, s = stacking, n = non-canonical, a = all, b = canonical + non-canonical)\n";
	cout << "\t (C = canonical or E = canonical + wobble (Extended), W = wobble, S = stacking, N = non-canonical, A = all, B = canonical + non-canonical)\n";
	cout << "\t you can choose more than one interaction for example -c abcns";
	cout << "\n";
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-v or --1D: comparison with 1D array algorithm\n";
	cout << "\t for determining the number of edges you wants to be present in calculation of similarity scores";
	cout << "\tyou can choose one of the following optinos (1, 2, 3, or 5 edges and faces)\n";
	cout << "\t\t Cans     --> 1, 3, or 5\n";
	cout << "\t\t NonCans  --> 3 or 5\n";
	cout << "\t\t Stacks   --> 2 0r 5\n";
	cout << "\t\t Wobbles  --> 1, 3, or 5\n";
	cout << "\t\t BasePairs--> 3 or 5\n";
	cout << "\t\t All      --> 5\n";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";

	cout << "Usage for clarna output files: \n\t./2D_N -r <referenceFile> -q <queryFile> -p <pdbFile> -c <requested_interaction> --1D <number_involved_edges> -o [outputName]\n";
	cout << "Example: \n\t./2D_N  -r samples/ClaRNARef.out -q samples/ClaRNAQuery.out -p samples/sample.pdb -c ebwn -o results/sample.csv --1D 3\n";
	exit(EXIT_FAILURE);
}

void printRequestedInteractions(string requestedInteractions)
{
	auto itA = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'A' || ch == 'a');});
	if(itA != cend(requestedInteractions))
	{
		cout << "\tAll" << endl;
	}
	
	auto itC = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'C' || ch == 'c');});
	auto itE = std::find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'E' or ch == 'e');});
	if(itC != cend(requestedInteractions) || itE != cend(requestedInteractions))
	{
		if(itE != cend(requestedInteractions))
		{
			cout << "\tCanonical + wobble" << endl;
		}
		else
		{
			cout << "\tCanonical" << endl;
		}
	}
	
	auto itN = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'N' || ch == 'n');});
	if(itN != cend(requestedInteractions))
	{
		cout << "\tNon-canonical" << endl;
	}
	
	auto itW = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'W' || ch == 'w');});
	if(itW != cend(requestedInteractions))
	{
		cout << "\tWobbles" << endl;
	}
	
	auto itS = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'S' || ch == 's');});
	if(itS != cend(requestedInteractions))
	{
		cout << "\tStacking" << endl;
	}
	
	auto itB = find_if(cbegin(requestedInteractions), cend(requestedInteractions), [](char ch){return (ch == 'B' || ch == 'b');});
	if(itB != cend(requestedInteractions))
	{
		cout << "\tAllBasePairs" << endl;
	}
}


int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
		exit(EXIT_FAILURE);
	}
    string refPath, queryPath, pdbFile, outFilename;
    int number_involved__edges_faces; 
    bool ClaRNAFlag, outFilenameFlag, refFlag, queryFlag, sequenceFlag, pdbFlag, vectorFlag;
    ClaRNAFlag = outFilenameFlag = refFlag = queryFlag = sequenceFlag  =  pdbFlag = vectorFlag = false;

    const struct option long_opts[] = {
    	{"query", required_argument, nullptr,'q'},
    	{"reference", required_argument, nullptr, 'r' },
    	{"pdb", required_argument, nullptr, 'p'},
    	{"clarna_out", required_argument, nullptr, 'c'},
    	{"1D", required_argument, nullptr, 'v'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	string requestedInteractions = "acnsw";
	int optind = 0;
	
	
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore_2D_N\n";
	cout << "===============================\n";

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"q:r:p:c:o:h", long_opts, &optind);
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
			case 'p':
				pdbFlag = true;
				pdbFile=optarg; 
				cout << "Input PDB file: " << pdbFile << endl;
				break;
			case 'r':
				refFlag = true;
				refPath = optarg;
				cout << "Reference path: " << optarg << endl;
				break;
			case 'c':
				ClaRNAFlag = true;
				requestedInteractions = optarg;
				cout << "Comparison of ClaRNA results will be done\n";
				cout << "Your requested interactions: \n";
				printRequestedInteractions(requestedInteractions);
				break;
			case 'v':
				vectorFlag = true;
				number_involved__edges_faces = stoi(optarg);
				cout << "number of involved edges and faces are: " << number_involved__edges_faces << endl;
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
		exit(EXIT_FAILURE);
	}
	
	if(!queryFlag)
	{
		cout << "Query file is required\n";
		showHelp();
		exit(EXIT_FAILURE);
	}
	
	if(!vectorFlag)
	{
		number_involved__edges_faces = 5;
		cout << "number of involved endges and faces are: " << number_involved__edges_faces << endl;
	}
	
	if(!outFilenameFlag)
	{
		outFilename = "FinalResuts.csv";
	}
	
	if(ClaRNAFlag)
	{
		if(!pdbFlag)
		{
			cout << "PDB file for -c option is required\n";
			showHelp();
			exit(EXIT_FAILURE);
		}
		
		CLARNA::CompareStructures cs(requestedInteractions);
		
		//cs.set_isWobble_canonical(wobble_canonicalFlag);
		cs.set_number_involved_faces_edges(number_involved__edges_faces);
		
		vector<ConfusionMatrixTuple> vcmtClaRNA = cs.calcConfusionMatrix(refPath, queryPath, pdbFile);
		cout << "sequence: " << cs.findInteractionRef.seqWithSeparateChains << endl;
		cs.writeScores(vcmtClaRNA, outFilename);
	}
	else
	{
		cout << "Please add opton -c in your command line" << endl;
		cout << "This optoin is required" << endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	
	cout << "Done ;)" << endl;

	
	
	return 0;
}
