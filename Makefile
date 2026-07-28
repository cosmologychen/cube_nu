#OPTIONS+=-CB
OPTIONS+=-DPID
OPTIONS+=-DZIPX
#OPTIONS+=-DHALOFIND
#OPTIONS+=-DSPEEDTEST

MODFILE:=$(wildcard *.f90)
OBJFILE:= Green.o  power_nu.o $(addprefix ,$(notdir $(MODFILE:.f90=.o))) cicpower_global.o

all: main.x
	@echo "done"
main.x: $(OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(OBJFILE) -o $@ $(FFTFLAG)
power_test.x: $(OBJFILE)
	@echo "Link files:"
	$(FC) $(XFLAG) $(OPTIONS) $(OBJFILE) -o $@ $(FFTFLAG)

parameters.o: Makefile basic_functions.f08
variables.o: parameters.o
pencil_fft.o: parameters.o
pencil_fft_global.o: parameters.o
cicpower_global.o: ./neutrinos/cicpower_global.f90 variables.o cubefft.o pencil_fft.o pencil_fft_global.o parameters.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
power_nu.o: ./neutrinos/power_nu.f90 pencil_fft.o cicpower_global.o
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)
timestep.0: power_nu.o variables.o 
particle_mesh.o: variables.o pencil_fft.o cubefft.o power_nu.o
initialize.o: variables.o Green.o cubefft.o pencil_fft.o power_nu.o z_checkpoint.txt z_halofind.txt
finalize.o: variables.o cubefft.o pencil_fft.o power_nu.o
kick.o: power_nu.o variables.o cubefft.o pencil_fft.o

main.o: $(OBJFILE)
$(OBJFILE): variables.o

power_test.o: $(OBJFILE)
$(OBJFILE): variables.o

#parameters.o: ../parameters.f90
#	$(FC) $(OFLAG) $(OPTIONS) $<
Green.o: ./Green/Green.f90
	$(FC) $(OFLAG) $(OPTIONS) $<
%.o: %.f90 Makefile
	$(FC) $(OFLAG) $(OPTIONS) $< -o $@ $(FFTFLAG)

clean:
	rm -vf *.mod *.o *.out *.err *.x *.tar hostfile_* *~
