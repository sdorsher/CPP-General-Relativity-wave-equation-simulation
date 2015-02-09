ITNT = -I../../tnt_126
IJAMA = -I../../jama125

dg1D : main.o ReferenceElement.o
	g++ -lm -std=c++11 $(ITNT) $(IJAMA)  main.o ReferenceElement.o -o dg1D

main.o: main.cpp ReferenceElement.h TNT2.h
	g++ -lm -std=c++11 $(ITNT) $(IJAMA) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp ReferenceElement.h TNT2.h
	g++ -lm -std=c++11 $(ITNT) $(IJAMA) -c ReferenceElement.cpp

