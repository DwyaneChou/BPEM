include ./configure.BPEM

OBJS = BPEM.o \
       constant.o \
       control.o \
       initialization.o \
       solve.o \
       integration.o \
       subroutines.o \
       output.o \
       tools.o

all: EXE

EXE:  $(OBJS)
	$(F90) -o $(EXENAME) $(OBJS) \
	-L$(NETCDF)/lib -lnetcdf -lnetcdff $(OPT) $(FCFLAGS) $(CFLAGS)

control.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) control.f90

constant.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) constant.f90

tools.o :
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) tools.f90

subroutines.o : control.o constant.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) subroutines.f90

integration.o : control.o constant.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) integration.f90

initialization.o : control.o tools.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) initialization.f90

solve.o : control.o subroutines.o integration.o output.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) solve.f90

output.o : control.o tools.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) output.f90

BPEM.o : control.o subroutines.o integration.o initialization.o solve.o output.o tools.o
	$(F90) -c $(OPT) $(FCFLAGS) $(CFLAGS) -L$(NETCDF)/lib -lnetcdf -lnetcdff $(INCLUDE) BPEM.f90

clean:
	$(RM) *.o *.mod *.exe *.a
