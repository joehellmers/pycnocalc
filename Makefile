F90COMP = gfortran
OBJS = globalvars.o logging.o vectors.o nucleon_interactions.o folding_potential.o astrophysics.o configuration.o utilities.o general_nuclear.o integration.o system_functions.o
MODULES = globalvars.mod logging.mod vectors.mod nucleon_interactions.mod folding_potential.mod astrophysics.mod configuration.mod utilities.mod general_nuclear.mod integration.mod system_functions.mod
SWITCHES = -ffree-line-length-none -O3
OUTDIR = /home/hellmersjl/bin

pycnocalc: $(MODULES)
	$(F90COMP) pycnocalc.f90 $(OBJS) -o $(OUTDIR)/pycnocalc $(SWITCHES)
	cp pycnocalc.cfg $(OUTDIR)
	chmod ugo-x $(OUTDIR)/pycnocalc.cfg

folding_potential.mod: constants.mod nucleon_interactions.mod general_nuclear.mod 
	$(F90COMP) -c folding_potential.f90 $(SWITCHES)


nucleon_interactions.mod: constants.mod vectors.mod 
	$(F90COMP) -c nucleon_interactions.f90 $(SWITCHES)


globalvars.mod: constants.mod
	$(F90COMP) -c globalvars.f90 $(SWITCHES)

vectors.mod: constants.mod
	$(F90COMP) -c vectors.f90 $(SWITCHES)

constants.mod: 
	$(F90COMP) -c constants.f90 $(SWITCHES)

logging.mod: constants.mod
	$(F90COMP) -c logging.f90 $(SWITCHES)

astrophysics.mod: constants.mod globalvars.mod
	$(F90COMP) -c astrophysics.f90 $(SWITCHES)

utilities.mod: constants.mod
	$(F90COMP) -c utilities.f90 $(SWITCHES)

configuration.mod: constants.mod logging.mod utilities.mod
	$(F90COMP) -c configuration.f90 $(SWITCHES)

general_nuclear.mod: constants.mod
	$(F90COMP) -c general_nuclear.f90 $(SWITCHES)

system_functions.mod: constants.mod astrophysics.mod logging.mod folding_potential.mod configuration.mod utilities.mod general_nuclear.mod integration.mod
	$(F90COMP) -c system_functions.f90 $(SWITCHES)

integration.mod: constants.mod
	$(F90COMP) -c integration.f90 $(SWITCHES)


clean: 
	rm *.mod
	rm *.o
	rm $(OUTDIR)/pycnocalc


