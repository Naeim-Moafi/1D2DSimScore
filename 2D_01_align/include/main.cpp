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
	cout << "1D2DSimScore1D_01_align\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-i --inputFile \n";
	//cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	cout << "-a or --align: for comparison structure inside the input file\n";
	cout << "\t (the type of the interaction must be chosen)\n";
	cout << "\t (c = canonical or e = canonical + wobble (extended), n = non-canonical, a = all)\n";
	cout << "\t (C = canonical or E = canonical + wobble (extended), N = non-canonical, A = all)\n";
	cout << "\t you can choose more than one interaction for example -d cna";
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	

	cout << "Usage for blast hit files:\n\t./2D_01_align -i <inputFile> -a <requested_interactions> -o [outputName]\n";
	cout << "Example:\n\t./2D_01_align -i samples/blast_example.txt -a can -o results/outputTest.csv\n";
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
		exit(EXIT_FAILURE);
	}
	
    string inpPath, outFilename; 
    bool  outFilenameFlag, inpFileFlag, wobble_canonicalFlag, blastHitFlag;
    outFilenameFlag = inpFileFlag = wobble_canonicalFlag = blastHitFlag = false;

    const struct option long_opts[] = {
    	{"inputFile", required_argument, nullptr,'i'},
    	{"align", no_argument,  nullptr, 'a'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
    
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore1D_01_align\n";
	cout << "===============================\n";
	string requestedInteractions = "acn";
	int optind = 0;

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"i:a:o:wh", long_opts, &optind);
        if (opt == -1)
        {
            break;
		}
        
        switch (opt)
        {
			case 'i':
				inpFileFlag = true;
				inpPath = optarg;
				cout << "Input path: " << optarg << endl;
				break;
				
			case 'a':
				blastHitFlag = true;
				requestedInteractions = optarg;
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
    if(!inpFileFlag)
    { 
		cout << "input file is required\n";
		showHelp();
	}
	
	
	if(!outFilenameFlag)
	{
		outFilename = "FinalResuts.csv";
	}
	
	if(blastHitFlag)
	{
		BLAST::CompareStructures cs;
		cs.requestedInteractions = requestedInteractions;
		//cs.findInteractions.set_isWobble_canonical(wobble_canonicalFlag);
		//cs.findInteractionRef.set_withSeq(sequenceFlag);
		//cs.findInteractionQuery.set_withSeq(sequenceFlag);
		
		auto vcmt = cs.calcConfusionMatrix(inpPath);
		cs.writeScores(vcmt, outFilename);
		
	}
	else
	{
		cout << "Please add opton -a in your command line" << endl;
		cout << "This optoin is obligatory for integrity of usage whole software" << endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	cout << "Done ;)" << endl;

	
	
	return 0;
}
