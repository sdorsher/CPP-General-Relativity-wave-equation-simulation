CXX = g++
#CXX = icpc
IGEN = -I/Users/sdorsher/Documents/Diener -I/home/sdorsher
ITNT = -I/home/sdorsher/tnt -I/Users/sdorsher/Documents/Diener/tnt -I/home/knarf/codes/dorsher/libtnt/libtnt-1.2.6/src/ -I /home/knarf/codes/dorsher/libjama/libjama-1.2.4/src
#IGEN = -I/home/sdorsher
#ITNT = -I/home/sdorsher/tnt
LCONF = -L/Users/sdorsher/utils/lib/ -lconfig++


dg1D : main.o ReferenceElement.o Grid.o Evolution.o globals.o ConfigParams.o DiffEq.o CharacteristicFlux.o Modes.o HyperboloidalCoords.o
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) main.o  ReferenceElement.o Grid.o Evolution.o globals.o ConfigParams.o DiffEq.o CharacteristicFlux.o Modes.o HyperboloidalCoords.o $(LCONF) -o dg1D
	uname | grep -q Linux || install_name_tool -change /usr/local/lib/libconfig++.9.dylib $(HOME)/utils/lib/libconfig++.9.dylib dg1D

main.o: main.cpp GridFunction.h GridFunction.tpp ReferenceElement.h VectorGridFunction.h VectorGridFunction.tpp TwoDVectorGridFunction.h TwoDVectorGridFunction.tpp Evolution.h DiffEq.h TNT2.h ConfigParams.h Modes.h HyperboloidalCoords.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) $(LCONF) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp ReferenceElement.h globals.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c ReferenceElement.cpp

#GridFunction.o: GridFunction.cpp GridFunction.h TNT2.h
#	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c GridFunction.cpp

#VectorGridFunction.o: VectorGridFunction.cpp VectorGridFunction.h GridFunction.h TNT2.h
#	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c VectorGridFunction.cpp

Grid.o: Grid.cpp Grid.h ReferenceElement.h GridFunction.h GridFunction.tpp TNT2.h DiffEq.h CharacteristicFlux.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c Grid.cpp

Evolution.o: Evolution.cpp Evolution.h GridFunction.h GridFunction.tpp TwoDVectorGridFunction.h TwoDVectorGridFunction.tpp ReferenceElement.h TNT2.h ConfigParams.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) $(LCONF) -c Evolution.cpp

globals.o: globals.cpp globals.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c globals.cpp

ConfigParams.o: ConfigParams.cpp ConfigParams.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) $(LCONF) -c ConfigParams.cpp

DiffEq.o: DiffEq.cpp DiffEq.h ConfigParams.h Grid.h CharacteristicFlux.h VectorGridFunction.h VectorGridFunction.tpp Modes.h HyperboloidalCoords.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) $(LCONF) -c DiffEq.cpp

CharacteristicFlux.o: CharacteristicFlux.cpp CharacteristicFlux.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c CharacteristicFlux.cpp

HyperboloidalCoords.o: HyperboloidalCoords.cpp HyperboloidalCoords.h globals.h 
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c HyperboloidalCoords.cpp

Modes.o: Modes.cpp Modes.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c Modes.cpp

clean:
	rm -f *.o dg1D
