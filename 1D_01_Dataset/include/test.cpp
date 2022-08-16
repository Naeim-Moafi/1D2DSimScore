#include "CompareStructures.h"
using namespace std;
using namespace Binary_all;

int main(int argc, char* argv[])
{
	string refPath = argv[1];
	string scores = argv[2];
	string outputPath = argv[3];
	CompareStructures cs;
	cs.init_binaries(refPath);
	cs.m_requestedScores = scores;
	cs.writeScores(outputPath);
	
	return EXIT_SUCCESS;
}
