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
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	

	cout << "Usage:\n\t./2D_01_CMO -r <referenceFile> -q <queryFile> -o [outputName]";
	cout << "\nexample:\n\t./2D_01_CMO  -r samples/ref.map -q samples/query.map -o results/sampleTest.csv";
	cout << endl;
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
	cout << "1D2DSimScore_2D_01_CMO\n";
	cout << "===============================\n";
	
    string refPath, queryPath, seqFile, outFilename; 
    [[maybe_unused]]bool dotBracketFlag, ClaRNAFlag, outFilenameFlag, refFlag, queryFlag, sequenceFlag, vectorFlag, matrixFlag;
    dotBracketFlag = ClaRNAFlag = outFilenameFlag = refFlag = queryFlag = sequenceFlag = vectorFlag = matrixFlag = false;

    const struct option long_opts[] = {
    	{"query", required_argument, nullptr,'q'},
    	{"reference", required_argument, nullptr, 'r' },
    	{"outputName", required_argument, nullptr, 'o'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	string requestedInteractions = "acn";
	int optind = 0;

    while (true)
    {
    	const auto opt = getopt_long(argc, argv,"q:r:o:h", long_opts, &optind);
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
			case 'r':
				refFlag = true;
				refPath = optarg;
				cout << "Reference path: " << optarg << endl;
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
	CMO::CompareStructures cs;
	
	cs.readFiles(refPath, queryPath);
	auto cmt = cs.calcConfusionMatrix();
	cs.writeScores(cmt, outFilename);
	
	auto end = high_resolution_clock::now();
	auto duration_s = duration_cast<seconds>(end - start);
	auto duration_ms = duration_cast<milliseconds>(end - start);
	
	if(duration_s.count() == 0) cout << "running time: " << duration_ms.count() << " ms\n";
	else cout << "running time: " << duration_s.count() << " s\n";
	
	cout << "Done ;)" << endl;

	
	
	return 0;
}
