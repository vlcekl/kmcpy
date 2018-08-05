#FC = mpif90
FC = gfortran
#FFLAGS = -fopenmp -O2 -ffast-math -w
FFLAGS = -O2 -ffast-math -w
#FFLAGS = -O0 -ffast-math -w

FILES = rndmod.f90 lgnvt.f90

all: clean

lgnvt: $(FILES)
	$(FC) $(FFLAGS) -o $@ $(FILES)

clean: lgnvt
	rm *.mod
