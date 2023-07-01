 
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




#include <chrono>
#include "Structures.h"
#include <getopt.h>
#include <chrono>

using namespace std;
using namespace std::chrono;


void showHelp(){
	cout << endl;
	cout << "Required arguments:\n";
	cout << "Either separate file for each biomolecule with at least two of the following options:\n";
	cout << "\t-R --rna: for RNA pdb file\n";
	cout << "\t-D --dna: for DNA pdb file\n";
	cout << "\t-P --protein: for Protein pdb file\n";
	cout << "or all in one with:\n";
	cout << "\t-p or --pdb: for input pdb file (all in one)\n";
	cout << "\t and\n";
	cout << "\t-c or --complex: for type of the complex (eg. RP for (RNA-Protein complex) or RD(RNA-DNA hybride) or RPD(for all possible interactions))";
	cout << "-t or --calculation_type: to determine the type of the calculatoin\n";
	cout << "\t INTRA --> intra-chain (contact between residues in the same chain)\n";
	cout << "\t INTER --> inter-chain (contact between residues in different chains)\n";
	cout << "\t ALL   --> all contacts\n";
	//cout << "-o or --outputName : name of output file (the format of the output is csv, so it is better for output file have the csv extention name.csv)\n";
	cout << endl;
	cout << "Optional arguments:\n";
	cout << "-d or --distance: the program considers wobble as canonical interactions (default is 6.0 A)\n";
	cout << "-n or --number_atoms_in_contact: it indicates least number of atoms in contact of each residue to consider two residue are in a interaction\n";
	cout << "-o or --outputName: name of output file (path + basename, and program will add its extension to that)\n";
	cout << "-m or --map: to create and visualize residue-residue contact maps\n";
	cout << "-S or --selection: to select only part of the strcuters for calculation\n";
	cout << "-C or --contact-checker: to check only specific contacts (for one chain structures only)\n";

	cout << "Usage:\n\t./ContactExtractor -p samples/3wbm.pdb -c RP -d 7 -n 3 -t INTER -m -o results/3wbm\n";
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	if(argc < 4){
		showHelp();
	}
	
    [[maybe_unused]]string rnaPath, dnaPath, proteinPath, pdbPath, outFilename, requested_complex, calculation_type, selection_path;
	[[maybe_unused]] string ci_path; 
    double distance = 6.0;
    int number_atoms_in_contact = 2;
    bool  outFilenameFlag, rnaFlag, dnaFlag, proteinFlag, pdbFlag, distanceFlag, numberFlag, complexFlag, mapFlag, typeFlag, selection_flag;
	bool ci_flag = false;
     outFilenameFlag = rnaFlag = dnaFlag = proteinFlag = pdbFlag = distanceFlag = numberFlag = complexFlag = mapFlag = selection_flag = false;

    const struct option long_opts[] = {
    	{"rna", required_argument, nullptr,'R'},
    	{"dna", required_argument, nullptr,'D'},
    	{"protein", required_argument, nullptr,'P'},
    	{"pdb", required_argument,  nullptr, 'p'},
    	{"distance", required_argument, nullptr, 'd'},
    	{"complex", required_argument, nullptr, 'c'},
    	{"number_atoms_in_contact", required_argument, nullptr, 'n'},
    	{"outputName", required_argument, nullptr, 'o'},
    	{"calculation_type", required_argument, nullptr, 't'},
    	{"selction", required_argument, nullptr, 's'},
    	{"map", no_argument, nullptr, 'm'},
		{"--contact-checker", required_argument, nullptr, 'C'},
    	{"help", no_argument, nullptr, 'h'},
    	{0,0,0,0}
    };
    
	int optind = 0;

    while (true)
    {
        //const auto opt = getopt_long(argc, argv,"t:f:c::r::dmo", long_opts, nullptr);
    	const auto opt = getopt_long(argc, argv,"R:D:P:p:d:n:c:o:t:C:S:mh", long_opts, &optind);
        if (opt == -1)
        {
            break;
		}
        
        switch (opt)
        {
			case 'R':
				rnaFlag = true;
				rnaPath = optarg;
				requested_complex += "R";
				cout << "RNA path: " << optarg << endl;
				break;
			
			case 'D':
				dnaFlag = true;
				dnaPath = optarg;
				cout << "DNA path: " << optarg << endl;
				requested_complex += "D";
				break;
			
			case 'P':
				proteinFlag = true;
				proteinPath = optarg;
				cout << "Protein path: " << optarg << endl;
				requested_complex += "P";
				break;
			case 'p':
				pdbFlag = true;
				pdbPath = optarg;
				cout << "input pdb path: " << optarg << endl;
				break;
			case 'c':
				complexFlag = true;
				requested_complex = optarg;
				break;
			
			case 'd':
				distanceFlag = true;
				distance = atof(optarg);
				break;
			
			case 'n':
				numberFlag = true;
				number_atoms_in_contact = atoi(optarg);
				break;
			
			case 'o':
				outFilenameFlag = true;
				outFilename = optarg;
				break;
			
			case 't':
				typeFlag = true;
				calculation_type = optarg;
				break;
			
			case 'm':
				mapFlag = true;
				break;
				
			case 'S':
				selection_flag = true;
				selection_path = optarg;
				break;
			case 'C': 
				ci_flag = true;
				ci_path = optarg;
				break;

			case 'h':
			case '?': 
			default:
				showHelp();
				break;
		}
    }
    
    // cheking input files
    if(!rnaFlag && !rnaFlag && !dnaFlag && !proteinFlag && !pdbFlag)
    { 
		cout << "you missed all required arguments\n";
		showHelp();
	}

	if(pdbFlag && !complexFlag)
	{
		cout << "you missed -c or --complex which is required\n";
		showHelp();
	}

	if((rnaFlag && complexFlag) || (dnaFlag && complexFlag) || (proteinFlag && complexFlag))
	{
		cout << "you can not use -c with -R, -D, -P at the same time\n";
		showHelp();
	}
	

	
	if(!outFilenameFlag)
	{
		outFilename = "resuts";
	}

	Structures st;
	if(typeFlag)
	{
		st.set_calculation_type(calculation_type);
	}
	else
	{
		cout << "Calculation type is required\n";
		showHelp();
	}
	
	if(selection_flag)
	{
		st.update = true;
		st.selection_path = selection_path;
	}
	st.set_is_plot_requested(mapFlag);

	st.set_distance_threshold(distance);
	st.set_number_atoms_in_contact(number_atoms_in_contact);

	auto start = chrono::high_resolution_clock::now();
	if(pdbFlag)
	{
		st.set_requested_molecules(requested_complex);
		st.readPDB_allStructures(pdbPath);
		//st.calc_distance_contact();
		//st.write_binary_format(outFilename);
	}
	else
	{
		auto itR = std::find_if(requested_complex.begin(), requested_complex.end(), [](char c){return (c == 'R' || c == 'r');});
		auto itD = std::find_if(requested_complex.begin(), requested_complex.end(), [](char c){return (c == 'D' || c == 'd');});
		auto itP = std::find_if(requested_complex.begin(), requested_complex.end(), [](char c){return (c == 'P' || c == 'p');});

		if(itR != requested_complex.end())
		{
			st.readPDB_rna(rnaPath);
		}

		if(itD != requested_complex.end())
		{
			st.readPDB_dna(dnaPath);
		}

		if(itP != requested_complex.end())
		{
			st.readPDB_protein(rnaPath);
		}

		//st.create_all_seq();
		//st.calc_distance_contact();
		//st.write_binary_format(outFilename);
	}
	
	//st.read_restraint_file("restraint.txt");
	//getchar();

	st.calc_distance_contact();
	st.write_binary_format(outFilename);
	st.write_as_dot_bracket(outFilename);
	if(mapFlag) st.write_map(outFilename);

	if(ci_flag)
	{
		st.check_contacts_interest(ci_path, outFilename + ".csv");
	}


	auto end = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(end - start);
	
	
	cout << "it took " << duration.count() << "s" << endl;
	
	
	cout << "Done ;)" << endl;	
	return 0;
}
