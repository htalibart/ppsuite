TARGET_LIB=ppalign_solver.so

#Sources files

CSRC = graph_apurva.cpp dp_mat_apurva.cpp problem_apurva.cpp bandb.cpp Domains.cpp main_ppalign.cpp
COBJ = $(CSRC:.cpp=.o)

#C compilator
CC = g++ #-fopenmp #-pg
#flags
#commented line is for gcc >= 4.3 
#OPTIM = -O3 -march=native -pipe -fomit-frame-pointer -ffast-math -mfpmath=sse,387
#OPTIM = -O3 -march=core2 -pipe -fomit-frame-pointer -ffast-math -mfpmath=sse,387
OPTIM = -O3 -march=core2 -ffast-math #-DNUM_THREADS=1 -fopenmp
DEBUG = -O0 -Wall -g -Wextra -ffast-math #-DNUM_THREADS=1 -fopenmp -pg
#CXXFLAGS = $(DEBUG)
CXXFLAGS = $(OPTIM) -fPIC -O2 -g
#CXXFLAGS = $(OPTIM) -fPIC -Wall -Wextra -O2 -g
LFLAGS = -lm -shared

#main work
$(TARGET_LIB): $(COBJ)
	$(CC) $(STDFLAG) -o $@ $(COBJ) $(LFLAGS)
	@echo
	@echo "    *******************************"
	@echo "    * Lib Lagrangian_Branch&Bound *"
	@echo "    *******************************"
	@echo "    -> Compilation complete !!!"	
	@echo
	
Bounds: main_boundcomp.o
	$(CC) -o $@ main_boundcomp.o $(LFLAGS)

mrf: graph_apurva.o dp_mat_apurva.o problem_apurva.o bandb.o
	$(CC) -I msgpack-c/mrf.h
	$(CC) -o $@ graph_apurva.o dp_mat_apurva.o problem_apurva.o bandb.o main_ppalign.o $(LFLAGS)


#remove .o
clean:
	rm -f *.o *~ *.aux *.log *.dvi *.toc

#
cleanall:
	rm -f *.o *~ *.aux *.log *.dvi *.toc A_purva documentation.pdf

	
doc: 
	doxygen doxygen.conf

doc_pdf:
#rules for all .o
graph_apurva.o : main.h
dp_mat_apurva.o : main.h
problem.o : main.h
problem_apurva.o : main.h
bandb.o : main.h
main.o: main.h
