srcdir = src
objdir = .obj
debug ?= 0 

FC = gfortran
FLAGS        = -I$(objdir) -J$(objdir)
LFLAGS		 = 

ifeq ($(debug), 1)
    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
else
    DFLAGS   = -O3
endif


## Individual libraries or modules ##
$(objdir)/greb.original.model.o: $(srcdir)/greb.original.model.f90
	@mkdir -p $(@D)
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/greb.o: $(srcdir)/greb.f90
	@mkdir -p $(@D)
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## greb executable
greb-original: $(objdir)/greb.original.model.o
	mkdir -p output
	$(FC) $(DFLAGS) $(FLAGS) -o greb-original $^ $(srcdir)/greb.original.shell.web-public.f90 $(LFLAGS)

greb: $(objdir)/greb.o
	mkdir -p output
	$(FC) $(DFLAGS) $(FLAGS) -o greb $^ $(LFLAGS)


clean:
	rm -f greb greb-original $(objdir)/*.o $(objdir)/*.mod
	rm -rf $(objdir)
