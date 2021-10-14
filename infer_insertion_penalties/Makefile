TARGET_LIB=inference_lib.so


CSRC = main_infer_insertion_penalties.cpp
COBJ = $(CSRC:.cpp=.o)

CC=g++
OPTIM = -O3 -march=core2 -ffast-math
CXXFLAGS = $(OPTIM) -fPIC -O2 -g
LFLAGS = -lm -shared


$(TARGET_LIB): $(COBJ)
	$(CC) $(STDFLAG) -o $@ $(COBJ) $(LFLAGS)

clean:
	rm -f *.o *~ *.aux *.log *.dvi *.toc


main_infer_insertion_penalties.o: main_infer_insertion_penalties.h
