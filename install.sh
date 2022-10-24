#!/bin/bash
# A script for installing different modules of 1D2DSimScore

install_1D_01() {
	echo "start installing 1D_01 package...";
	cd 1D_01;
	make > /dev/null ;
	cp "1D_01" ../bin
	cd ..;
	echo "1D_01 package is installed";
}

install_1D_01_Dataset() {
	echo "start installing 1D_01_Dataset package...";
	cd 1D_01_Dataset;
	make > /dev/null;
	cp "1D_01_Dataset" ../bin
	cd ..;
	echo "1D_01_Dataset package is installed";
}

install_2D_01() {
	echo "start installing 2D_01 package...";
	cd 2D_01;
	make > /dev/null ;
	cp "2D_01" ../bin
	cd ..;
	echo "2D_01 is installed";
}


install_2D_01_Dataset() {
	echo "start installing 2D_01_Dataset package...";
	cd 2D_01_Dataset;
	make > /dev/null;
	cp "2D_01_Dataset" ../bin
	cd ..;
	echo "2D_01_Dataset package is installed";
}

install_2D_01_align() {
	echo "start installing 2D_01_align package...";
	cd 2D_01_align;
	make > /dev/null ;
	cp "2D_01_align" ../bin
	cd ..;
	echo "2D_01_align package is installed";
}

install_2D_01_CMO() {
	echo "start installing 2D_01_CMO package...";
	cd 2D_01_CMO;
	make > /dev/null ;
	cp "2D_01_CMO" ../bin
	cd ..;
	echo "2D_01_CMO package is installed";
}

install_2D_N() {
	echo "start installing 2D_N package...";
	cd 2D_N;
	make > /dev/null ;
	cp "2D_N" ../bin
	cd ..;
	echo "2D_N package is installed";
}

install_2D_N_Dataset() {
	echo "start installing 2D_N_Dataset package...";
	cd 2D_N_Dataset;
	make > /dev/null;
	cp "2D_N_Dataset" ../bin
	cd ..;
	echo "2D_N_Dataset package is installed";
}









install_all() {
	install_1D_01;
	install_1D_01_Dataset;
	install_2D_01;
	install_2D_01_Dataset;
	install_2D_01_align;
	install_2D_01_CMO;
	install_2D_N;
	install_2D_N_Dataset;
}

clean_all(){
	cd 1D_01;
	make clean > /dev/null;
	cd ..;
	cd 1D_01_Dataset;
	make clean > /dev/null;
	cd ..;
	cd 2D_01_align;
	make clean > /dev/null;
	cd ..
	cd 2D_01_CMO;
	make clean > /dev/null;
	cd ..;
	cd 2D_N;
	make clean > /dev/null;
	cd ..;
	cd 2D_N_Dataset;
	make clean > /dev/null;
	cd ..;
	cd 2D_01;
	make clean > /dev/null;
	cd ..;
	cd 2D_01_Dataset;
	make clean > /dev/null;
	cd ..;
	rm -rf bin;
}

showHelp(){
	echo "For installing you can write name of the packages as arguments separated by spaces"
	echo "or All for installing all packages";
	echo
	echo "For unistallling you can choose clean";
}


mkdir -p bin;
if [ "$#" -lt 1 ];
then
	showHelp;
	exit 0;
fi

for arg in $@
do 
	case $arg in
		1D_01) install_1D_01;;
		1D_01_Dataset) install_1D_01_Dataset;;
		2D_01) install_2D_01;;
		2D_01_Dataset) install_2D_01_Dataset;;
		2D_01_align) install_2D_01_align;;
		2D_01_CMO) install_2D_01_CMO;;
		2D_N) install_2D_N;;
		2D_N_Dataset) install_2D_N_Dataset;;
		all)install_all;;
		clean)clean_all;;
		*) echo "$arg is an unknown package"; showHelp; exit 0;;
	esac
done
