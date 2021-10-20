TARGET_LIB=cpp_rescale_lib.so


CSRC = cpp_rescale.cpp
COBJ = $(CSRC:.cpp=.o)

CC=g++
OPTIM = -O3 -march=core2 -ffast-math
CXXFLAGS = $(OPTIM) -fPIC -O2 -g
LFLAGS = -lm -shared


$(TARGET_LIB): $(COBJ)
	$(CC) $(STDFLAG) -o $@ $(COBJ) $(LFLAGS)

clean:
	rm -f *.o *~ *.aux *.log *.dvi *.toc


cpp_rescale.o: cpp_rescale.h
