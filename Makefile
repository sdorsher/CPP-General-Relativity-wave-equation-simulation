CXX = g++
#CXX = icpc
IGEN = -I/Users/sdorsher/Documents/Diener -I/home/sdorsher
ITNT = -I/home/sdorsher/tnt -I/Users/sdorsher/Documents/Diener/tnt -I/home/knarf/codes/dorsher/libtnt/libtnt-1.2.6/src/ -I /home/knarf/codes/dorsher/libjama/libjama-1.2.4/src
#IGEN = -I/home/sdorsher
#ITNT = -I/home/sdorsher/tnt


dg1D : main.o ReferenceElement.o
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) main.o ReferenceElement.o -o dg1D

main.o: main.cpp ReferenceElement.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp ReferenceElement.h TNT2.h
	$(CXX) -g -lm -std=c++11 $(ITNT) $(IGEN) -c ReferenceElement.cpp
