#CXX = g++
CXX = icpc
IGEN = -I/home/sdorsher
ITNT = -I/home/sdorsher/tnt

dg1D : main.o ReferenceElement.o
	$(CXX) -lm -std=c++11 $(ITNT) $(IGEN) main.o ReferenceElement.o -o dg1D

main.o: main.cpp ReferenceElement.h TNT2.h
	$(CXX) -lm -std=c++11 $(ITNT) $(IGEN) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp ReferenceElement.h TNT2.h
	$(CXX) -lm -std=c++11 $(ITNT) $(IGEN) -c ReferenceElement.cpp
