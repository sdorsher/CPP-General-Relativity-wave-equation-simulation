CXX = g++
#CXX = icpc
IGEN = -I/Users/sdorsher/Documents/Diener -I/home/sdorsher
ITNT = -I/home/sdorsher/tnt -I/Users/sdorsher/Documents/Diener/tnt -I/home/knarf/codes/dorsher/libtnt/libtnt-1.2.6/src/ -I /home/knarf/codes/dorsher/libjama/libjama-1.2.4/src
#IGEN = -I/home/sdorsher
#ITNT = -I/home/sdorsher/tnt


dg1D : main.o GridFunction.o ReferenceElement.o VectorGridFunction.o Grid.o Evolution.o
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) main.o GridFunction.o ReferenceElement.o VectorGridFunction.o Grid.o Evolution.o globals.o -o dg1D

main.o: main.cpp GridFunction.h ReferenceElement.h VectorGridFunction.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp ReferenceElement.h globals.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c ReferenceElement.cpp

GridFunction.o: GridFunction.cpp GridFunction.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c GridFunction.cpp

VectorGridFunction.o: VectorGridFunction.cpp VectorGridFunction.h GridFunction.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c VectorGridFunction.cpp

Grid.o: Grid.cpp Grid.h ReferenceElement.h GridFunction.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c Grid.cpp

Evolution.o: Evolution.cpp Evolution.h GridFunction.h VectorGridFunction.h ReferenceElement.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c Evolution.cpp

globals.o: globals.cpp globals.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c globals.cpp