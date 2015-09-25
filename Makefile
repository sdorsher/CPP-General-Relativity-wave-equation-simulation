CXX = gcc
#CXX = icpc

#Spine
IGEN = -I/home/sdorsher -I/home/sdorsher/libconfig-1.5/install/include
ITNT = -I/home/sdorsher/tnt -I/home/sdorsher/jama
LCONF = -L/home/sdorsher/libconfig-1.5/execinstall/lib/ -lconfig++
ESRC = /home/sdorsher/scalar1deffectivesource
LGSL = `pkg-config --libs gsl`


#Both
LCPP = -lstdc++
FLGS = -g -lm -std=c++11 -O3 -p


#Steven's Mac
#IGEN = -I/Users/sdorsher/Documents/Diener -I/home/sdorsher
#ITNT = -I/home/sdorsher/tnt -I/Users/sdorsher/Documents/Diener/tnt -I/home/knarf/codes/dorsher/libtnt/libtnt-1.2.6/src/ -I /home/knarf/codes/dorsher/libjama/libjama-1.2.4/src
#LCONF = -L/Users/sdorsher/utils/lib/ -lconfig++
#ESRC = /Users/sdorsher/Documents/Diener/Scalar1DEffectiveSource/scalar1deffectivesource
#LGSL=?


#LGSL = -L/usr/local/packages/gsl-1.16-intel/lib -lgsl -lgslcblas



#dg1D : main.o ReferenceElement.o Grid.o Evolution.o globals.o ConfigParams.o DiffEq.o CharacteristicFlux.o Modes.o HyperboloidalCoords.o  orbit.o $(ESRC)/EffectiveSource.o $(ESRC)/WignerDMatrix.o namespaces.o numerics.o source_interface.o
dg1D : main.o ReferenceElement.o Grid.o Evolution.o globals.o ConfigParams.o DiffEq.o CharacteristicFlux.o Modes.o HyperboloidalCoords.o  orbit.o namespaces.o numerics.o source_interface.o WriteFile.o
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) main.o  ReferenceElement.o Grid.o Evolution.o globals.o ConfigParams.o orbit.o DiffEq.o CharacteristicFlux.o Modes.o HyperboloidalCoords.o namespaces.o  $(ESRC)/EffectiveSource.o $(ESRC)/WignerDMatrix.o numerics.o source_interface.o WriteFile.o $(LCONF) $(LGSL) $(LCPP) -o dg1D
	uname | grep -q Linux || install_name_tool -change /usr/local/lib/libconfig++.9.dylib $(HOME)/utils/lib/libconfig++.9.dylib dg1D

main.o: main.cpp GridFunction.h GridFunction.tpp ReferenceElement.h VectorGridFunction.h VectorGridFunction.tpp TwoDVectorGridFunction.h TwoDVectorGridFunction.tpp Evolution.h DiffEq.h TNT2.h ConfigParams.h Modes.h HyperboloidalCoords.h namespaces.h orbit.h source_interface.h numerics.h WriteFile.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) $(LCONF) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp ReferenceElement.h globals.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -I$(ESRC) -c ReferenceElement.cpp

#GridFunction.o: GridFunction.cpp GridFunction.h TNT2.h
#	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c GridFunction.cpp

#VectorGridFunction.o: VectorGridFunction.cpp VectorGridFunction.h GridFunction.h TNT2.h
#	$(CXX) $(FLGS) $(ITNT) $(IGEN) -c VectorGridFunction.cpp

Grid.o: Grid.cpp Grid.h ReferenceElement.h GridFunction.h GridFunction.tpp TNT2.h DiffEq.h CharacteristicFlux.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c Grid.cpp

Evolution.o: Evolution.cpp Evolution.h GridFunction.h GridFunction.tpp TwoDVectorGridFunction.h TwoDVectorGridFunction.tpp ReferenceElement.h TNT2.h ConfigParams.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) $(LCONF) -c Evolution.cpp

globals.o: globals.cpp globals.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c globals.cpp

ConfigParams.o: ConfigParams.cpp ConfigParams.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) $(LCONF) -c ConfigParams.cpp

DiffEq.o: DiffEq.cpp DiffEq.h ConfigParams.h Grid.h CharacteristicFlux.h VectorGridFunction.h VectorGridFunction.tpp Modes.h HyperboloidalCoords.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) $(LCONF) -c DiffEq.cpp

CharacteristicFlux.o: CharacteristicFlux.cpp CharacteristicFlux.h TNT2.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c CharacteristicFlux.cpp

HyperboloidalCoords.o: HyperboloidalCoords.cpp HyperboloidalCoords.h globals.h 
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c HyperboloidalCoords.cpp



Modes.o: Modes.cpp Modes.h TwoDVectorGridFunction.h ConfigParams.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c Modes.cpp

namespaces.o: namespaces.cpp namespaces.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c namespaces.cpp

orbit.o: orbit.cpp orbit.h ConfigParams.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) $(LCONF) -c orbit.cpp

numerics.o: numerics.cpp numerics.h 
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c numerics.cpp

source_interface.o: source_interface.cpp source_interface.h Modes.h numerics.h namespaces.h Grid.h VectorGridFunction.h GridFunction.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c source_interface.cpp

WriteFile.o: WriteFile.cpp WriteFile.h TwoDVectorGridFunction.h Grid.h DiffEq.h
	$(CXX) $(FLGS) $(ITNT) $(IGEN) -I$(ESRC) -c WriteFile.cpp

clean:
	rm -f *.o dg1D
