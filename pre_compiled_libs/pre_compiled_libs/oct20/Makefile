include make.inc.TOKS_IPP

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

LLIB=tglf_lib

ifdef FANN_ROOT
   # neural-net file
   NN=tglf_nn_TM.fann
else
   # stub file
   NN=tglf_nn_TM
endif

OBJECTS = tglf_isnan.o      \
          tglf_modules.o    \
          tglf_pkg.o        \
          tglf_allocate.o   \
          tglf_deallocate.o \
          tglf_startup.o    \
          tglf_hermite.o    \
          tglf_inout.o      \
          tglf_setup_geometry.o \
          tglf_LS.o         \
          tglf_eigensolver.o \
          tglf_geometry.o   \
          tglf_matrix.o     \
          tglf_max.o        \
          tglf_interface.o  \
          tglf_error.o      \
          tglf_isinf.o      \
          tglf_shutdown.o   \
          tglf_read_input.o \
          tglf_multiscale_spectrum.o \
          tglf_kygrid.o     \
          tglf_run.o        \
          tglf_TM.o         \
          tglf_init_mpi.o   \
          tglf_run_mpi.o    \
          tglf_TM_mpi.o     \
          $(NN).o

.SUFFIXES : .o .f90 .F90

all: $(LLIB).a 

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

tglf_run_mpi.o: tglf_run.F90
	$(FC) -o tglf_run_mpi.o $(FMATH) $(FFLAGS) -DMPI_TGLF=1 -c $<

.f90.o:
	$(FC) $(FMATH) $(FFLAGS) -qopennp $(NN_LIB) -c $<
.F90.o:
	$(FC) $(FMATH) $(FFLAGS) $(NN_LIB) -c $<

clean:
	rm -f *.o *.lst $(LLIB).a 
	rm -f ${GACODE_ROOT}/modules/tglf*.mod

-include make.inc.TOKS_IPP
