CXXFLAGS =	-O2 -Wall -pipe -march=native -std=c++20

CXX = g++-11

LIBS = -pthread

OBJS =   include/Point3D.o include/ReadPDB.o include/Structures.o include/main.o #  include/test.o

TARGET = include/Point3D.o	include/ReadPDB.o include/Structures.o
  

all: contactExtractor

contactExtractor: $(OBJS)
	    $(CXX) -o ContactExtractor $(OBJS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET) ContactExtractor
