#FCOMPFLAGS    =	-O0 -g -C -fpe0 -D nb_num -traceback
#FCOMPFLAGS    =	-O2  -g -C -fpe0 -traceback -check noarg_temp_created 
# FCOMPFLAGS    =	-Os -march=native 
FCOMPFLAGS    =	-Os
#FCOMPFLAGS    =	 -pg -tpp7 -O3  -align  -axN -tpp7 
#FCOMPFLAGS    =	 -pg  -O2  
# CCOMPFLAGS    =	-Os -march=native
CCOMPFLAGS    =	-Os

FFLAGS        =	$(FCOMPFLAGS)
CFLAGS        = $(CCOMPFLAGS)
LDFLAGS       =	$(FCOMPFLAGS)

LD            =	gfortran
FC            = gfortran

MAKEFILE      =	Makefile
PROGRAM       =	scattcg 



FSOURCE = $(wildcard src/*.F)
OBJS1 = $(patsubst %.F, %.o, $(FSOURCE))
OBJS2 = $(patsubst src/%, build/%, $(OBJS1))

INCS = $(wildcard src/*.inc)

build/%.o: src/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 

$(PROGRAM): clean $(OBJS2) $(INCS)
	$(LD) $(LDFLAGS) $(OBJS2) -o $(PROGRAM)  

clean:
	rm -f $(OBJS2) $(PROGRAM)

all: $(PROGRAM)
