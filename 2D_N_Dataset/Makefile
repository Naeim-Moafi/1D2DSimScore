CXXFLAGS =	-O2 -Wall -pipe -march=native -std=c++20

CXX = g++-11

OBJS   = include/FindInteractions.o include/lib/FindInteractions.o include/lib/ReadPDB.o include/CompareStructures.o include/lib/CompareStructures.o include/lib/SimilarityScores.o include/main.o 

TARGET = include/FindInteractions.o include/lib/FindInteractions.o include/lib/ReadPDB.o include/CompareStructures.o include/lib/CompareStructures.o include/lib/SimilarityScores.o
  

all: CompareStructures

CompareStructures: $(OBJS)
	    $(CXX) -o 2D_N_Dataset $(OBJS) 

clean:
	rm -f $(OBJS) $(TARGET) 2D_N_Dataset
