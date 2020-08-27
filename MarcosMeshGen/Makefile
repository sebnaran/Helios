CXXFLAGS = -O

PRJ = _05

2D : 
	$(CXX) $(CXXFLAGS) ./SRC-MAIN/main$(PRJ).cc $(INC) -c
	$(CXX) main$(PRJ).o -o main$(PRJ)
	rm main$(PRJ).o

clean:
	make clean_SRC
	make clean_EXE
	make clean_MSH

clean_SRC:
	mv SRC*/*~ TRASH/.

clean_EXE:
	mv main_* TRASH/.

clean_MSH:
	rm ./MESH-PS/*.ps
