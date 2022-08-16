Author:	
		Seyed Naeim Moafinejad: snmoafinejad@iimcb.gov.pl
# ContactMapExtractor
# 1. Functionality
Generating distance contact map

# 2. Installation
## Prerequisites
---
- g++ v11 (g++-11)
- gnuplot
---

For installation use make command.
For unistall the software you can use make clean command. 

# 3. Usage

Required arguments:

Either separate file for each biomolecule with at least two of the following options:

	-R --rna: for RNA pdb file
	-D --dna: for DNA pdb file
	-P --protein: for Protein pdb file

or all in one with:

	-p or --pdb: for input pdb file (all in one)
	 and
	-c or --complex: for type of the complex (eg. RP for (RNA-Protein complex) or RD(RNA-DNA hybride) or RPD(for all possible interactions))-t or --calculation_type: to determine the type of the calculatoin
	 INTRA --> intra-chain (contact between residues in the same chain)
	 INTER --> inter-chain (contact between residues in different chains)
	 ALL   --> all contacts

Optional arguments:

-d or --distance: the program considers wobble as canonical interactions (default is 6.0 A)

-n or --number_atoms_in_contact: it indicates least number of atoms in contact of each residue to consider two residue are in a interaction

-o or --outputName: name of output file (path + basename, and program will add its extension to that)

-m or --map: to create and visualize residue-residue contact maps

Usage:
	
	./ContactExtractor -p samples/3wbm.pdb -c RP -d 7 -n 3 -t INTER -m -o results/3wbm

