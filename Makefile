mods = mod_read_tools.f90 mod_random.f90 mod_netdata.f90
program = dynamics.f90
comp = gfortran
c?=
ifeq (${c},1)
	#ch=-check all -traceback # ifort
	ch=-g -fbacktrace -fbounds-check -ffpe-trap=zero,overflow,underflow -Wall -fcheck=all # gfortran
endif

.PHONY: dynamics clean
default: dynamics

dynamics: clean
	rm -f mod_read_tools.mod mod_random.mod mod_netdata.mod ifport.mod
	$(comp) $(mods) $(program) -o dynamics ${ch}
	
clean:
	rm -f mod_read_tools.mod mod_random.mod mod_netdata.mod ifport.mod