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
	cout << "1D2DSimScore_2D_01_Dataset\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-i --input (it could be a directory  with files for different structures or single file with all informations) \n";
	cout << "-D or --dot_bracket_all: for comparison structure inside the input file\n";
	cout << "-S or --scores: choose the score for calculation of similarity scores al vs all\n";
	cout << "\t MCC        --> Mathew Correlation Coeficient\n";
	cout << "\t FSCORE     --> FSCore\n";
	cout << "\t JINDEX     --> Jaccard Index\n";
	cout << "\t FMINDEX    --> FMIndex\n";
	cout << "\t PRECISION  --> Precision\n";
	cout << "\t RECALL     --> Recall\n";
	cout << "\t SPECIFICITY--> Specificity\n";
	cout << "\t BA         --> Balanced Accuracy\n";
	cout << "\t FOR        --> False Omission Rate\n";
	cout << "\t PT         --> Prevalence Threshold\n";
	cout << "\t CSI        --> Critical Success Index\n";
	cout << "\t MK         --> MarKedness\n";
	cout << "\t JBINDEX    --> Janusz Bujnicki Index\n";
	cout << "\t(you can ask for a score with upper or lower case or any combination of them)\n";
	
	cout << endl;
	cout << "Optional arguments:\n";
	//cout << "-w or --wobble_canonical : the program considers wobble as canonical interactions\n";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	cout << "-m or --2D: if you prefer to calculate the similarity scores with 2D algorithm\n";
	cout << "-e or --extension: extension of the files in your dataset, you can change the extension when you use a directory as your input (default: .SS)\n";
	//cout << "-f or --flags4scores: List all required scores. (Eg -f \"TP FSCORE GSCORE JINDEX\"). Default is FSCORE\n";
	//cout << "-m or --meta_data_write : To write meta data of TP, FP, FN to *.tpstr file\n";
	//cout << "-g or --setting : To set the ranges of residues to consider for comparision\n";
	

	cout << "Usage for 2D_01_Dataset with a file:\n\t./2D_01_Dataset -i <inputFile> -D -S <requested_scores_separated_with_comma> -o [outputName]\n";
	cout << "Exmple for 2D_01_Dataset with a file:\n\t./2D_01_Dataset -i samples/AllInOne.SS_all  -D -S MCC,Fscore,for,jInDeX --2D -o results/test.gsm\n";
	cout << "Usage for 2D_01_Dataset with a folder:\n\t./2D_01_Dataset -i <inputFolder> -D -S <requested_scores_separated_with_comma> -o [outputName]\n";
	cout << "Exmple for 2D_01_Dataset with a folder:\n\t./2D_01_Dataset -i samples/SSs  -D -S MCC,Fscore,for,jInDeX --2D -o results/test.gsm\n";
	cout << "Exmple for 2D_01_Dataset with a folder:\n\t./2D_01_Dataset -i samples/dbns  -D -S MCC,Fscore,for,jInDeX --2D -o results/test.gsm\n";
	cout << "Exmple for 2D_01_Dataset with a folder and different extension:\n\t./2D_01_Dataset -i samples/dbns  -D -S MCC,Fscore,for,jInDeX --2D -o results/test.gsm -e .dbn\n";
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
		exit(EXIT_FAILURE);
	}
	
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore_2D_01_Dataset\n";
	cout << "===============================\n";
	
    string inpPath, outFilename, extension = ".SS"; 
    bool  outFilenameFlag, inpFlag, scoreFlag, d_allFlag, matrixFlag, extensionFlag;
    outFilenameFlag = inpFlag = scoreFlag = d_allFlag = matrixFlag = extensionFlag = false;

    const struct option long_opts[] = {
    	{"input", required_argument, nullptr,'i'},
    	{"dot_bracket_all", no_argument,  nullptr, 'D'},
    	{"2D", no_argument,  nullptr, 'm'},
    	{"scores", required_argument, nullptr, 'S'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"extension", required_argument, nullptr, 'e'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	string requested_scores;// = "123456";
	int optind = 0;

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"i:S:o:e:Dmh", long_opts, &optind);
        if (opt == -1)
        {
            break;
		}
        
        switch (opt)
        {
			case 'i':
				inpFlag = true;
				inpPath = optarg;
				cout << "Input path: " << optarg << endl;
				break;
				
			case 'S':
				scoreFlag = true;
				requested_scores = optarg;
				break;
			
			case 'D':
				d_allFlag = true;
				break;
			
			case 'o':
				outFilenameFlag = true;
				outFilename = optarg;
				break;
			
			case 'e':
				extensionFlag = true;
				extension = optarg;
				break;
				
			case 'm':
				matrixFlag = true;
				break;
				
			case 'h': 
			case '?': 
			default:
				showHelp();
				break;
		}
    }
    
    // cheking input files
    if(!inpFlag)
    { 
		cout << "input file is required\n";
		showHelp();
	}
	
	if(!outFilenameFlag)
	{
		outFilename = "FinalResuts.gsm";
	}
	
	if(d_allFlag)
	{
		SS_all::CompareStructures cs;
		cs.findInteractions.set_extension(extension);
		cs.m_requestedScores = requested_scores;
		if(matrixFlag)
		{
			cs.set_is_matrix_prefered(true);
		}
		cs.readFile(inpPath);
		cs.writeScores(outFilename);
	}
	else
	{
		cout << "Please add opton -D in your command line" << endl;
		cout << "This optoin is obligatory for integrity of usage whole software" << endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	cout << "Done ;)" << endl;

	
	
	return 0;
}
