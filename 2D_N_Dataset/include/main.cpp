 
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
	cout << "1D2DSimScore_2D_N_Dataset\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-i or --inputs: reference file (it can be bracket dot notation of SS or ClaRNA.out depending on option )\n";
	cout << "-C or --clarna_all:  for comparison *.out files with each other (all vs all)\n";
	cout << "\t (the type of the interaction must be chosen (there is no defualt type))\n";
	cout << "\t (c = canonical or e = canonical + wobble (extended), w = wobble, s = stacking, n = non-canonical, a = all, b = canonical + non-canonical)\n";
	cout << "\t (C = canonical or E = canonical + wobble (Extended), W = wobble, S = stacking, N = non-canonical, A = all, B = canonical + non-canonical)\n";
	cout << "\t you can choose more than one interaction for example -c abcns";
	cout << "-S or --scores: choose the score for calculation of similarity scores al vs all\n";
	cout << "\t MCC        --> Mathew Correlation Coeficient\n";
	cout << "\t FSCORE     --> F1 Score\n";
	cout << "\t JINDEX     --> Jaccard Index\n";
	cout << "\t FMINDEX    --> FM Index\n";
	cout << "\t PRECISION  --> Precision\n";
	cout << "\t RECALL     --> Recall\n";
	cout << "\t SPECIFICITY--> Specificity\n";
	cout << "\t BA         --> Balanced Accuracy\n";
	cout << "\t FOR        --> False Omission Rate\n";
	cout << "\t PT         --> Prevalence Threshold\n";
	cout << "\t CSI        --> Critical Success Index\n";
	cout << "\t MK         --> MarKedness\n";
	cout << "\t JBINDEX    --> Bujnicki Index\n";
	cout << "\t(you can ask for a score with upper or lower case or any combination of them)\n";
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
	cout << "\t\t All      --> 3 or 5\n";
	cout << "-e or --extension : extension of the files in your dataset, you can change the extension when you use a directory as your input (default: .out)\n";
	cout << "-o or --outputName : name of output file (the format of the output is gsm, so it is better for output file have the csv extention name.gsm)\n";

	cout << "Usage for 2D_N_Dataset: \n\t./2D_N_Dataset -i <inputDirectory> -C <requested_interactions> --1D <number_of_involved_edges> -S <requested_scores_separated_with_comma> -o [outputName]\n";
	cout << "Examples for 2D_N_Dataset: \n\t./2D_N_Dataset -i samples/outs -C cansb --1D 5 -S fscore,mcc,mk,csi,jindex,recall -o results/sample_test.gsm\n";
	cout << "Examples for 2D_N_Dataset with different extension: \n\t./2D_N_Dataset -i samples/contacts -C cansb --1D 5 -S fscore,mcc,mk,csi,jindex,recall -o results/sample_test_contact.gsm -e .contact\n";
	
	exit(EXIT_FAILURE);
}


int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
		exit(EXIT_FAILURE);
	}
	
    string inputPath, outFilename, extension;
    int number_involved__edges_faces = 5; 
    bool ClaRNAFlag, outFilenameFlag, inputFlag, vectorFlag, scoreFlag, extFlag;
    ClaRNAFlag = outFilenameFlag = inputFlag = vectorFlag = scoreFlag = extFlag = false;

    const struct option long_opts[] = {
    	{"inputs", required_argument, nullptr,'i'},
    	{"clarna_all", required_argument, nullptr, 'C'},
    	{"1D", required_argument, nullptr, 'v'},
    	{"scores", required_argument, nullptr, 'S'},
		{"extension", required_argument, nullptr, 'e'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	string requestedInteractions = "acnsw";
	string requestedScores;
	int optind = 0;
	
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore_2D_N_Dataset\n";
	cout << "===============================\n";

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"i:C:S:o:e:wh", long_opts, &optind);
        if (-1 == opt)
        {
            break;
		}
        
        switch (opt)
        {
			case 'i':
				inputFlag = true;
				inputPath = optarg; 
				cout << "input directory: " << optarg << endl;
				break;
				
			case 'C':
				ClaRNAFlag = true;
				requestedInteractions = optarg;
				cout << "Comparison of ClaRNA results will be done (all vs all)\n";
				break;
				
			case 'S':
				scoreFlag = true;
				requestedScores = optarg;
				break;
				
			case 'v':
				vectorFlag = true;
				number_involved__edges_faces = stoi(optarg);
				cout << "number of involved edges and faces are: " << number_involved__edges_faces << endl;
				break;
			
			case 'e':
				extFlag = true;
				extension = optarg;
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
    if(!inputFlag)
    { 
		cout << "input directory is required\n";
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
		CLARNA_all::CompareStructures cs;
		if(extFlag) cs.findInteractions.set_extension(extension);
		if(!scoreFlag)
		{
			cout << "-S is required\n";
			showHelp();
			exit(EXIT_FAILURE);
		}
		cs.m_requestedScores = requestedScores;
		cs.requestedInteractions = requestedInteractions;
		//cs.set_isWobble_canonical(wobble_canonicalFlag);
		cs.set_number_involved_faces_edges(number_involved__edges_faces);
		cs.findInteractions.readInputFiles(inputPath);
		cout << "requested interactions and scores\n";
		cs.writeScores(outFilename);
	}
	else
	{
		cout << "-C is required" <<endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	
	std::filesystem::path path(outFilename);
	cout << "your results are in " << path.parent_path() << endl;
	cout << "Done ;)" << endl;

	
	
	return 0;
}
