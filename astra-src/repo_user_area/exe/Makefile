ifndef AWD
$(error setenv/export AWD first)
endif

# Paths

LIB_DIR = $(AWD)/lib
TSK_DIR = $(AWD)/bin
ASTRA_DIR = $(AWD)/src/astra
MOD_DIR = $(AWD)/mod
OBJ_DIR = $(AWD)/obj
TMP_DIR = $(AWD)/src/tmp
XPR_DIR = $(AWD)/xpr
EQU = $(AWD)/equ/$(MODEL)
EXP = $(AWD)/exp/$(FEXP)
LIBSBR   = $(LIB_DIR)/libsbr.a
LIBFNC   = $(LIB_DIR)/libfnc.a
LIBASTRA = $(LIB_DIR)/libastra.a
LIBMISC  = $(LIB_DIR)/libmisc.a
LIBNBI   = $(LIB_DIR)/libnbi.a
LIBFEQIS = $(LIB_DIR)/libfeqis.a
LIBGRAPH = $(LIB_DIR)/libgraph.a
MOD_MISC  = $(MOD_DIR)/misc
MOD_ASTRA = $(MOD_DIR)/astra
MOD_FEQIS = $(MOD_DIR)/feqis
MOD_GRAPH = $(MOD_DIR)/graph
MOD_NBI   = $(MOD_DIR)/nbi
AOBJ_DIR = $(OBJ_DIR)/astra
FOBJ_DIR = $(OBJ_DIR)/feqis
GOBJ_DIR = $(OBJ_DIR)/graph
MOBJ_DIR = $(OBJ_DIR)/misc
NOBJ_DIR = $(OBJ_DIR)/nbi

FC_FLAGS += -I$(AWD) # to incldue fml in sbr/*.f90
LD_FLAGS += -I$(AWD) # to incldue fml in src/tmp/*.f90

# parser-generated Fortran routines

TMP_SRC = \
$(TMP_DIR)/associate_pointers.f90 \
$(TMP_DIR)/inivar.f90 \
$(TMP_DIR)/detvar.f90 \
$(TMP_DIR)/postep.f90 \
$(TMP_DIR)/ininam.f90 \
$(TMP_DIR)/init_converge_step.f90 \
$(TMP_DIR)/eqns_inc.f90

# graphics yes/no

PLATFORM := $(shell $(AWD)/get_platform)

ifdef X11
# Debian-family multiarch puts libX11 in /usr/lib/<arch>-linux-gnu, RPM-family
# in /usr/lib64.  Probe for the library instead of naming a directory per
# architecture: this covers aarch64, which used to need its own branch, and
# x86_64 Debian, which had none and so failed to link the GUI build.
   MULTIARCH_LIBDIR = /usr/lib/$(shell uname -m)-linux-gnu
   ifeq ($(PLATFORM),darwin)
      X11_LIB = -L/opt/X11/lib -lX11
   else ifneq ($(wildcard $(MULTIARCH_LIBDIR)/libX11.so*),)
      X11_LIB = -L$(MULTIARCH_LIBDIR) -lX11
   else
      X11_LIB = -L/usr/lib64 -lX11
   endif
   X11_FLAG = -DX11
   GRAPH_FLAG = -lgraph
   TMP_SRC += $(TMP_DIR)/astra_out.f90 $(TMP_DIR)/set_graph_names.f90
else
   X11_LIB =
   X11_FLAG =
   GRAPH_FLAG =
endif

# ASTRA objects

AMODS = \
$(AOBJ_DIR)/json_vars.o \
$(AOBJ_DIR)/debugger.o \
$(AOBJ_DIR)/machine_config.o \
$(AOBJ_DIR)/scrunch2d.o \
$(AOBJ_DIR)/stella_module.o \
$(AOBJ_DIR)/vmec_indata.o \
$(AOBJ_DIR)/parse_utils.o \
$(AOBJ_DIR)/emeq.o \
$(AOBJ_DIR)/scalars.o \
$(AOBJ_DIR)/read_input.o \
$(AOBJ_DIR)/time_functions.o \
$(AOBJ_DIR)/status.o \
$(AOBJ_DIR)/transport_solver.o \
$(AOBJ_DIR)/standard_functions.o \
$(AOBJ_DIR)/cpu_usage.o \
$(AOBJ_DIR)/parameters_a2equil.o \
$(AOBJ_DIR)/surface_contours.o \
$(AOBJ_DIR)/flux_avg_imp.o \
$(AOBJ_DIR)/json_rw.o \
$(AOBJ_DIR)/gs_solver.o \
$(AOBJ_DIR)/a2vmec.o \
$(AOBJ_DIR)/a2dkes.o \
$(AOBJ_DIR)/metrics.o \
$(AOBJ_DIR)/auxiliary.o \
$(AOBJ_DIR)/set_x_data.o

AOBJS = $(AOBJ_DIR)/ipc_control.o
AINCLUDE = -I$(MOD_MISC) -I$(MOD_ASTRA) -I$(JSON_INC)

# FEQIS objects

FMODS = \
$(FOBJ_DIR)/transport2fbe.o \
$(FOBJ_DIR)/modules_core.o \
$(FOBJ_DIR)/modules_feqis.o \
$(FOBJ_DIR)/feqis_tools.o \
$(FOBJ_DIR)/fbe_core.o \
$(FOBJ_DIR)/pbe_core.o \
$(FOBJ_DIR)/circuit.o \
$(FOBJ_DIR)/feqis_solvers.o

FOBJS = \
$(FOBJ_DIR)/solver_inversion_matrix_gsef.o \
$(FOBJ_DIR)/feqis_main.o \
$(FOBJ_DIR)/PHI_EQ_2d_PBE.o \
$(FOBJ_DIR)/jacobians.o \
$(FOBJ_DIR)/build_2dgrid.o 

FINCLUDE = -I$(MOD_MISC) -I$(MOD_FEQIS) -I$(JSON_INC)

# graph objects

GMODS = $(GOBJ_DIR)/graph_utils.o $(GOBJ_DIR)/gui_interaction.o
GOBJS = $(GOBJ_DIR)/Astra2XW.o $(GOBJ_DIR)/dialogx.o
GINCLUDE = -I$(MOD_MISC) -I$(MOD_ASTRA) -I$(MOD_GRAPH)

# misc objects

MMODS = \
$(MOBJ_DIR)/char_manip.o \
$(MOBJ_DIR)/dbl2char.o \
$(MOBJ_DIR)/imas_ids.o \
$(MOBJ_DIR)/numerical_tools.o \
$(MOBJ_DIR)/pi_const.o

# NBI objects

NMODS = \
$(NOBJ_DIR)/cross_sections.o \
$(NOBJ_DIR)/nbstatus.o \
$(NOBJ_DIR)/nbicom.o \
$(NOBJ_DIR)/nbibce.o \
$(NOBJ_DIR)/nb_injection.o

NINCLUDE = -I$(MOD_MISC) -I$(MOD_ASTRA) -I$(MOD_NBI)

# Linking terms

#FOR STELLARATOR
#LIB = -lnbi -lsbr -lnbi -L/usr/lib64 -lX11 $(QLKNN_LIB) $(TGLF_LIB) $(TRAVIS_LIB) $(TORBEAM_LIB) $(RABBIT_LIB) $(GRAPH_FLAG) -lastra -lfnc -lfeqis -lmisc $(SPIDER_LIB) $(JSON_LIB) -L$(LIB_DIR) $(MKL_LIB) -lrt $(TORIC_LIB) -lstdc++
#LIB +=  -Wl,--start-group $(LIBSTELL_LIBDIR) -Wl,--end-group  -lrt  $(NETCDF_LIB) 

MKL_LIB = -Wl,--start-group $(MKL_LIBDIR)/libmkl_intel_lp64.a $(MKL_LIBDIR)/libmkl_core.a $(MKL_LIBDIR)/libmkl_sequential.a -Wl,--end-group
RT_LIB = -lrt
ifeq ($(PLATFORM),darwin)
   MKL_LIB = -L$(OPENBLAS_LIBDIR) -lopenblas
   RT_LIB =
endif

LIB = -lsbr -lnbi $(X11_LIB) $(QLKNN_LIB) $(TGLF_LIB) $(TRAVIS_LIB) $(TORBEAM_LIB) $(RABBIT_LIB) $(GRAPH_FLAG) -lastra -lfnc -lfeqis -lmisc $(SPIDER_LIB) $(JSON_LIB) -L$(LIB_DIR) $(MKL_LIB) $(RT_LIB) $(TORIC_LIB)

AENV = $(AWD)/exe/astra_rc

ASTRA_SRC = \
$(ASTRA_DIR)/astra_main.F90 $(WORKFLOW)

PY_TMP_GEN = $(AWD)/pyparse/parser_main.py -equ $(EQU) -exp $(EXP)

FMLS    := $(wildcard $(AWD)/fml/*)
PYFILES := $(wildcard $(AWD)/pyparse/*.py)

# Make Commands

.PHONY: all
all: directories modules $(ATASK)

.PHONY: directories
directories:
	@mkdir -p $(TSK_DIR)
	@mkdir -p $(TMP_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(AWD)/ncdf_out
	@mkdir -p $(AWD)/dat

.PHONY: clean
clean:
	rm -f *.mod
	make -f $(AWD)/exe/Makexpr clean
	rm -rf $(MOD_DIR)
	rm -rf $(TSK_DIR)
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)

.PHONY: modules
modules:
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(MMODS)"                       sub=misc  modules
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(AMODS)" "INC_FLAGS=$(AINCLUDE)" sub=astra modules
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(FMODS)" "INC_FLAGS=$(FINCLUDE)" sub=feqis modules
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(GMODS)" "INC_FLAGS=$(GINCLUDE)" sub=graph modules
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(NMODS)" "INC_FLAGS=$(NINCLUDE)" sub=nbi   modules

# Linking, create ASTRA executable

$(ATASK): $(LIBGRAPH) $(LIBASTRA) $(LIBMISC) $(LIBNBI) $(LIBFEQIS) $(SPIDER_LIB) $(LIBSBR) $(LIBFNC) $(AENV) $(EQU) $(ASTRA_SRC) $(PYFILES) $(FMLS)
	$(PY_TMP_GEN)
	$(FC) $(LD_FLAGS) -o $(ATASK) $(ASTRA_SRC) $(TMP_SRC) -I$(MOD_MISC) -I$(MOD_ASTRA) -I$(MOD_DIR)/sbr -I$(MOD_FEQIS) -I$(MOD_GRAPH) $(LIB) $(X11_FLAG)

# User defined sbr and fnc

$(LIBSBR): $(AWD)/sbr/*
	make -f $(AWD)/exe/Makesub sub=sbr

$(LIBFNC): $(AWD)/fnc/*.f90
	make -f $(AWD)/exe/Makesub sub=fnc

# Libs

$(LIBMISC): $(AWD)/src/misc/*
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(MMODS)" sub=misc lib

$(LIBASTRA): $(AWD)/src/astra/*
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(AMODS)" "OBJS=$(AOBJS)" "INC_FLAGS=$(AINCLUDE)" sub=astra lib

$(LIBNBI): $(AWD)/src/nbi/*
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(NMODS)" "INC_FLAGS=$(NINCLUDE)" sub=nbi lib

$(LIBFEQIS): $(AWD)/src/feqis/*
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(FMODS)" "OBJS=$(FOBJS)" "INC_FLAGS=$(FINCLUDE)" sub=feqis lib

$(LIBGRAPH): $(AWD)/src/graph/*
	make -f $(AWD)/exe/Makeastra "MODOBJS=$(GMODS)" "OBJS=$(GOBJS)" "INC_FLAGS=$(GINCLUDE)" sub=graph lib
