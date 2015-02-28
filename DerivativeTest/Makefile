CXX = icpc
IGEN = -I/home/sdorsher
ITNT = -I/home/sdorsher/tnt
IDG1D = -I/home/sdorsher/DG1D-CPP-development/DG1D-CPP-development

dconv : main.o ReferenceElement.o
	$(CXX) -lm -std=c++11 -g $(ITNT) $(IGEN) $(IDG1D) main.o ReferenceElement.o -o dconv

main.o: main.cpp ../ReferenceElement.h ../TNT2.h ../tnt_array1D_extn.h ../tnt_array2D_extn.h
	$(CXX) -lm -std=c++11 -g $(ITNT) $(IGEN) $(IDG1D) -c main.cpp

ReferenceElement.o: ../ReferenceElement.cpp ../ReferenceElement.h ../TNT2.h ../tnt_array1D_extn.h ../tnt_array2D_extn.h
	$(CXX) -lm -std=c++11 $(ITNT) $(IGEN) $(IDG1D) -c ../ReferenceElement.cpp
