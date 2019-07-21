F90COMP = gfortran
OBJS = globalvars.o \
	logging.o \
	vectors.o \
	nucleon_interactions.o \
	folding_potential.o \
	astrophysics.o \
	configuration.o \
	utilities.o \
	general_nuclear.o \
	mathutils.o \
	mathlinearalg.o \
	mathinterpolation.o \
	mathintegration.o \
	system_functions.o \
	cmdline.o \
	rate_calc.o \
	tests.o
MODULES = globalvars.mod \
	logging.mod \
	vectors.mod \
	nucleon_interactions.mod \
	folding_potential.mod \
	astrophysics.mod \
	configuration.mod \
	utilities.mod \
	general_nuclear.mod \
	mathutils.mod \
	mathlinearalg.mod \
	mathinterpolation.mod \
	mathintegration.mod \
	system_functions.mod \
	cmdline.mod \
	rate_calc.mod \
	tests.mod
SWITCHES = -ffree-line-length-none -O3 -fopenmp -fbounds-check
OUTDIR = ~/bin

pycnocalc: $(MODULES)
	$(F90COMP) pycnocalc.f90 $(OBJS) -o $(OUTDIR)/pycnocalc $(SWITCHES)
	cp pycnocalc.cfg $(OUTDIR)
	chmod ugo-x $(OUTDIR)/pycnocalc.cfg

folding_potential.mod: constants.mod nucleon_interactions.mod general_nuclear.mod configuration.mod 
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

configuration.mod: constants.mod logging.mod utilities.mod globalvars.mod
	$(F90COMP) -c configuration.f90 $(SWITCHES)

general_nuclear.mod: constants.mod mathintegration.mod
	$(F90COMP) -c general_nuclear.f90 $(SWITCHES)

system_functions.mod: constants.mod astrophysics.mod logging.mod folding_potential.mod configuration.mod utilities.mod general_nuclear.mod mathintegration.mod rate_calc.mod mathinterpolation.mod
	$(F90COMP) -c system_functions.f90 $(SWITCHES)

mathutils.mod: constants.mod
	$(F90COMP) -c mathutils.f90 $(SWITCHES)

mathlinearalg.mod: constants.mod
	$(F90COMP) -c mathlinearalg.f90 $(SWITCHES)

mathinterpolation.mod: constants.mod mathlinearalg.mod mathutils.mod
	$(F90COMP) -c mathinterpolation.f90 $(SWITCHES)

mathintegration.mod: constants.mod
	$(F90COMP) -c mathintegration.f90 $(SWITCHES)

cmdline.mod:
	$(F90COMP) -c cmdline.f90 $(SWITCHES)

rate_calc.mod: constants.mod folding_potential.mod astrophysics.mod mathintegration.mod logging.mod general_nuclear.mod
	$(F90COMP) -c rate_calc.f90 $(SWITCHES)

tests.mod: constants.mod rate_calc.mod general_nuclear.mod logging.mod mathlinearalg.mod mathinterpolation.mod
	$(F90COMP) -c tests.f90 $(SWITCHES)

dist:
	scp $(OUTDIR)/pycnocalc electra:~/bin

all: pycnocalc dist


clean:
	rm -rf *.mod
	rm -rf *.o
	rm -rf .DS_Store
	rm -rf $(OUTDIR)/pycnocalc
	rm -rf notebooks/.ipynb_checkpoints/
clean_local:
	rm *.mod
	rm *.o
