#CXX = icpc
CXX = g++
#IGEN = -I/Users/sdorsher/Documents/Diener
IGEN = -I/home/sdorsher
ITNT = -I/home/sdorsher/tnt
#ITNT = -I/Users/sdorsher/Documents/Diener/tnt
IDG1D = -I/home/sdorsher/DG1D-CPP-development/DG1D-CPP-development
#IDG1D = -I/Users/sdorsher/Documents/Diener/DG1D-CPP-development2/DG1D-CPP-development

dconv : main.o ReferenceElement.o
	$(CXX) -lm -std=c++11 -g $(ITNT) $(IGEN) $(IDG1D) main.o ReferenceElement.o -o dconv

main.o: main.cpp ../ReferenceElement.h ../TNT2.h ../tnt_array1D_extn.h ../tnt_array2D_extn.h
	$(CXX) -lm -std=c++11 -g $(ITNT) $(IGEN) $(IDG1D) -c main.cpp

ReferenceElement.o: ../ReferenceElement.cpp ../ReferenceElement.h ../TNT2.h ../tnt_array1D_extn.h ../tnt_array2D_extn.h
	$(CXX) -lm -std=c++11 $(ITNT) $(IGEN) $(IDG1D) -c ../ReferenceElement.cpp
