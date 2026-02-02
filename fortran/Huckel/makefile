F90 = gfortran
FFLAGS = -O2
LDFLAGS = -llapack 

OBJ = Huckel_mat.o subroutine.o diag.o

all: Huckel_mat

Huckel_mat: $(OBJ)
	$(F90) -o $@ $(OBJ) $(LDFLAGS)

%.o: %.f90
	$(F90) $(FFLAGS) -c $<

clean:
	rm -f *.o Huckel_mat

