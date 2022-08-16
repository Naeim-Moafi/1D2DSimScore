#include <algorithm>
#include "CompareStructures.h"
#include <getopt.h>

using namespace std;


void showHelp(){
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore1D_01\n";
	cout << "===============================\n";
	cout << "Required arguments:\n";
	cout << "-r or --reference: reference file\n";
	cout << "-q or --query: target file\n";
	//cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	cout << "-b or --binary: for comparison *.xo files with each other\n";
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	

	cout << "Usage:\n\t./1D_01 -r <referenceFile> -q <queryFile> -b -o [outputName]";
	cout << "\nexample:\n\t./1D_01  -r samples/ref.xo -q samples/query.xo -b -o results/sampleTest.csv";
	cout << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	if(argc < 2){
		showHelp();
	}
	
    string refPath, queryPath, outFilename; 
    bool binaryFlag, outFilenameFlag, refFlag, queryFlag;
    binaryFlag  = outFilenameFlag = refFlag = queryFlag = false;

    const struct option long_opts[] = {
    	{"query", required_argument, nullptr,'q'},
    	{"reference", required_argument, nullptr, 'r' },
    	{"binary", no_argument,  nullptr, 'b'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	string requestedInteractions = "acn";
	int optind = 0;
	
	cout << endl;
	cout << "===============================\n";
	cout << "1D2DSimScore1D_01\n";
	cout << "===============================\n";

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"q:r:o:bh", long_opts, &optind);
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
				
			case 'b':
				binaryFlag = true;
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
	
	if(binaryFlag)
	{
		Binary::CompareStructures cs;
		cs.readInputFiles(refPath, queryPath);
		ConfusionMatrixTuple cmt;
		
		cmt = cs.calcConfusionMatrix();
		cs.writeScores(cmt, outFilename);
	}
	else
	{
		cout << "Please add opton -b in your command line" << endl;
		cout << "This optoin is obligatory for integrity of usage of whole software" << endl;
		showHelp();
		exit(EXIT_FAILURE);
	}
	cout << "Done ;)" << endl;

	
	
	return 0;
}
