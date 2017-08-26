PROG =	SubChannelFlow

SRCS =	steam_nrs.f90 Main.f90 #ReactorCode.f90

OBJS =	steam_nrs.o Main.o #ReactorCode.o 

LIBS = -L/usr/local/lib -L/usr/lib 


CC = $(CC)
CFLAGS = -O0
#FC = $(FC)
#FFLAGS = -O0
FC      = gfortran
F90 = $(FC)
F90FLAGS =  -r8
F90FLAGS =  -g
LDFLAGS = -O0 -r8 #-cpp
F90FLAGS  = -ffast-math -frecord-marker=4 -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8 
LINKER  = g++   

all: $(PROG)

$(PROG):	$(OBJS)
		$(F90) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<
