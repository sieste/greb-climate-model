srcdir = src
objdir = .obj
debug ?= 0 

FC = gfortran
FLAGS        = -I$(objdir) -J$(objdir)
LFLAGS		 = 

ifeq ($(debug), 1)
    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
else
    DFLAGS   = -O
endif


## Individual libraries or modules ##
$(objdir)/greb.model.o: $(srcdir)/greb.model.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## greb executable
greb: $(objdir)/greb.model.o
	$(FC) $(DFLAGS) $(FLAGS) -o greb $^ $(srcdir)/greb.shell.web-public.f90 $(LFLAGS)

clean:
	rm -f greb $(objdir)/*.o $(objdir)/*.mod
