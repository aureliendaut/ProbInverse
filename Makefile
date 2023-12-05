# compiler
CC=gcc

# Optimization options
OPTIM_FLAG = -O3 -DNDEBUG -I  -std=c99 -Wall -lm
# Debug options
DEBUG_FLAG = -g -I -Wall -lm -std=c99
# Choose how to compile here
CXX_FLAGS = $(DEBUG_FLAG)

# The name of the executable file
PROG = chp 

# Sources to compile
SRC = main.c GC.c matrice.c function.c parametres.c

# The compilation line
$(PROG) : $(SRC)
	$(CC) $(SRC) $(CXX_FLAGS) -o $(PROG) 

# Compile and run
run : $(PROG)
	./$(PROG)

# Delete compilation files
clean :
	rm *.o *∼ *.dat *.log *.png $(PROG)

images : 
	gnuplot "gnuplot.txt"
	
images_exactes :
	gnuplot "gnuplot_exacte.txt"
