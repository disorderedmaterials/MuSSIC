#FCOMPFLAGS    =	-O0 -g -C -fpe0 -D nb_num -traceback
#FCOMPFLAGS    =	-O2  -g -C -fpe0 -traceback -check noarg_temp_created 
# FCOMPFLAGS    =	-Os -march=native 
FCOMPFLAGS    =	-Os -fbounds-check
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



FSOURCE = $(wildcard src_cg/*.F)
OBJS1 = $(patsubst %.F, %.o, $(FSOURCE))
OBJS2 = $(patsubst src_cg/%, build_cg/%, $(OBJS1))

INCS = $(wildcard src_cg/*.inc)

build_cg/%.o: src_cg/%.F
	$(FC) $(FFLAGS) -c $< -o $@ 

$(PROGRAM): clean $(OBJS2) $(INCS)
	$(LD) $(LDFLAGS) $(OBJS2) -o $(PROGRAM)  

clean:
	rm -f $(OBJS2) $(PROGRAM)

all: $(PROGRAM)
