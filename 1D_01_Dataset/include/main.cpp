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
	cout << "1D2DSimScore1D_01_Dataset\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-i or --input: input file of input directory\n";
	cout << "-q or --query: target file\n";
	cout << "-B or --binary_all: for comparison *.xo files with each other\n";
	cout << "-S or --scores: choose the score for calculation of similarity scores al vs all\n";
	cout << "\t MCC        --> Mathew Correlation Coeficient\n";
	cout << "\t FSCORE     --> F1-Score\n";
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
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-e or --extension : extension of the files in you dataset, your can change the extension when you use a directory as your input (default: .xo)";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extension name.csv)\n";
	

	cout << "Usage for 1D-01 dataset comparsion:\n\t./1D_01_Dataset -i <inputFile> -B -S <requested_scores_separated_with_comma> -o [outputName]\n";
	cout << "Exmple for 1D-01 dataset comparsion with a file:\n\t./1D_01_Dataset -i samples/sample1.xo  -B -S MCC,Fscore,for,jInDeX -o results/test.gsm\n";
	cout << "Usage for 1D-01 dataset comparsion with a folder:\n\t./1D_01_Dataset -i <inputFolder> -B -S <requested_scores_separated_with_comma> -o [outputName] -e [.extension]\n";
	cout << "Exmple for 1D-01 dataset comparsion with a folder:\n\t./1D_01_Dataset -i samples/XOs_dir  -B -S MCC,Fscore,for,jInDeX, JBINDEX -o results/test.gsm\n";
	cout << "Exmple for 1D-01 dataset comparsion with a folder and new extension:\n\t./1D_01_Dataset -i samples/01s_dir  -B -S MCC,Fscore,for,jInDeX, JBINDEX -o results/test.gsm -e .bin\n";
	cout << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
		exit(EXIT_FAILURE);
	}
	
    string inpPath, outFilename, extension; 
    bool binaryFlag, outFilenameFlag, inpFlag, scoreFlag, extFlag;
    binaryFlag = outFilenameFlag = inpFlag = scoreFlag = extFlag = false;

	const struct option long_opts[] = {
		{"input", required_argument, nullptr,'i'},
		{"binary_all", no_argument,  nullptr, 'B'},
		{"scores", required_argument, nullptr, 'S'},
		{"extension", required_argument, nullptr, 'e'},
		{"outputName", required_argument, nullptr, 'o'},
		{"help", no_argument, nullptr, 'h'},
		{0,0,0,0}
	};
    
	string requested_scores = "MCC";
	int optind = 0;
	
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore1D_01_Dataset\n";
	cout << "===============================\n";

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"i:S:o:e:Bh", long_opts, &optind);
        if (-1 == opt)
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
			
			case 'B':
				binaryFlag = true;
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
    if(!inpFlag)
    { 
		cout << "input file or input directory is required\n";
		showHelp();
	}
	
	
	
	if(!outFilenameFlag)
	{
		outFilename = "FinalResuts.gsm";
	}
	
	if(binaryFlag)
	{
		Binary_all::CompareStructures cs;
		cs.m_requestedScores = requested_scores;
		cs.init_binaries(inpPath);
		if(extFlag) cs.findInteractions.set_extension(extension);
		//cs.findInteractionRef.set_withSeq(sequenceFlag);
		//cs.findInteractionQuery.set_withSeq(sequenceFlag);
		
		cs.writeScores(outFilename);
	}
	else
	{
		cout << "Please add opton -d in your command line" << endl;
		cout << "This optoin is obligatory for integrity of usage whole software" << endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	cout << "Done ;)" << endl;

	
	
	return 0;
}
