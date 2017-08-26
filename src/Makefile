SRCS =	steam_nrs.f90 neutronics_fluids_subroutines.f90

OBJS =	steam_nrs.o neutronics_fluids_subroutines.o

#LIBS = -L/usr/local/lib -L/usr/lib 


F90 = gfortran
F90FLAGS  = -ffast-math -frecord-marker=4 -ffree-line-length-none -ffixed-line-length-none  -O3 -fdefault-real-8 
  

#all: $(PROG)

all: coupled fluids_only

coupled:	coupled_neutronics_fluids.o $(OBJS)
		$(F90) -o $@ coupled_neutronics_fluids.o $(OBJS) #$(LIBS)

fluids_only:	fluids_only.o $(OBJS)
		$(F90) -o $@ fluids_only.o $(OBJS) #$(LIBS)

#%.o:	%.f90
#	$(F90) -c $(F90FLAGS) $<

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

clean:
	rm -f coupled fluids_only $(OBJS) *.mod
