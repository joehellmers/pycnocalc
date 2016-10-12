F90COMP = gfortran
OBJS = globalvars.o logging.o vectors.o nucleon_interactions.o folding_potential.o
MODULES = globalvars.mod logging.mod vectors.mod nucleon_interactions.mod folding_potential.mod
SWITCHES = -ffree-line-length-none -O3
OUTDIR = /home/hellmersjl/bin

pycnocalc: $(MODULES)
	$(F90COMP) pycnocalc.f90 $(OBJS) -o $(OUTDIR)/pycnocalc $(SWITCHES)

folding_potential.mod: constants.mod nucleon_interactions.mod 
	$(F90COMP) -c folding_potential.f90 $(SWITCHES)


nucleon_interactions.mod: constants.mod vectors.mod 
	$(F90COMP) -c nucleon_interactions.f90 $(SWITCHES)


globalvars.mod: constants.mod globalvars.f90
	$(F90COMP) -c globalvars.f90 $(SWITCHES)

vectors.mod: vectors.f90
	$(F90COMP) -c vectors.f90 $(SWITCHES)

constants.mod: constants.f90
	$(F90COMP) -c constants.f90 $(SWITCHES)

logging.mod: logging.f90
	$(F90COMP) -c logging.f90 $(SWITCHES)

clean: 
	rm *.mod
	rm *.o
	rm $(OUTDIR)/pycnocalc


