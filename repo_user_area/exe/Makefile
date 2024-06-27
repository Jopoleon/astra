ifndef AWD
$(error setenv/export AWD first)
endif

ifndef WORKFLOW
WORKFLOW = $(FOR_DIR)/stepup.f90
endif

LIB_DIR = $(AWD)/lib
TSK_DIR = $(AWD)/bin
FOR_DIR = $(AWD)/main
MOD_DIR = $(AWD)/mod/for
OBJ_DIR = $(AWD)/obj/for
TMP_DIR = $(AWD)/tmp
NBI_INC = $(AWD)/mod/nbi
ATASK = $(TSK_DIR)/$(MODEL).exe
EQU  = $(AWD)/equ/$(MODEL)
LIBSBR = $(LIB_DIR)/libsbr.a
LIBFNC = $(LIB_DIR)/libfnc.a
LIBFOR = $(LIB_DIR)/libfor.a
LIBNBI = $(LIB_DIR)/libnbi.a
LIBFEQIS  = $(LIB_DIR)/libfeqis.a
LIBSPIDER = $(LIB_DIR)/libspider.a 

MKL_LIB = -Wl,--start-group $(MKL_LIBDIR)/libmkl_intel_lp64.a $(MKL_LIBDIR)/libmkl_core.a $(MKL_LIBDIR)/libmkl_sequential.a -Wl,--end-group
LIB = -lnbi -lsbr -lnbi -L/usr/lib64 -lX11 $(QLKNN_LIB) $(TGLF_LIB) $(TORB_LIB) $(RABBIT_LIB) -lfor -lfnc -lfeqis -lspider -L$(LIB_DIR) $(MKL_LIB) -lrt

AENV = $(AWD)/exe/astra_rc

AFC_FLAGS = $(CHECKALLCOMP) -qopenmp -g -traceback -o   

FOR_SRC = \
$(FOR_DIR)/astra_main.f90 \
$(WORKFLOW)

TMP_SRC = \
$(TMP_DIR)/astra_out.f90 \
$(TMP_DIR)/inivar.f90 \
$(TMP_DIR)/setvar.f90 \
$(TMP_DIR)/detvar_init.f90 \
$(TMP_DIR)/detvar.f90 \
$(TMP_DIR)/postep.f90 \
$(TMP_DIR)/subproc.f90 \
$(TMP_DIR)/ininam.f90 \
$(TMP_DIR)/init_converge_step.f90 \
$(TMP_DIR)/eqns_inc.f90

PY_TMP_GEN = $(AWD)/pyparse/parser_main.py -equ $(EQU)

FMLS    := $(wildcard $(AWD)/fml/*)
PYFILES := $(wildcard $(AWD)/pyparse/*.py)

# Make Commands
.PHONY: all
all: directories modules $(ATASK)

.PHONY: directories
directories:
	@mkdir -p $(TSK_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(AWD)/.res
	@mkdir -p $(AWD)/ncdf_out
	@mkdir -p $(AWD)/dat

.PHONY: clean
clean:
	rm -f *.mod
	rm -f $(MOD_DIR)/*.mod
	rm -f $(TSK_DIR)/*
	rm -rf $(OBJ_DIR)
	rm -rf $(LIB_DIR)

.PHONY: modules
modules:
	make -f $(AWD)/exe/Makefor modules
	make -f $(AWD)/exe/Makefeqis modules

# Linking, create ASTRA executable

$(ATASK): $(LIBFOR) $(LIBNBI) $(LIBFEQIS) $(LIBSPIDER) $(LIBSBR) $(LIBFNC) $(AENV) $(EQU) $(FOR_SRC) $(PYFILES) $(FMLS)
	$(PY_TMP_GEN)
	$(AFC) $(AFC_FLAGS) $(ATASK) $(FOR_SRC) $(TMP_SRC) -I$(MOD_DIR) -I$(NBI_INC) $(LIB)

# User defined sbr and fnc

$(LIBSBR): $(AWD)/sbr/*
	make -f $(AWD)/exe/Makesub sub=sbr

$(LIBFNC): $(AWD)/fnc/*.f90
	make -f $(AWD)/exe/Makesub sub=fnc

# Libs

$(LIBFOR): $(AWD)/src/for/*
	make -f $(AWD)/exe/Makefor lib

$(LIBNBI): $(AWD)/src/nbi/*
	cd $(AWD)/src; make -f $(AWD)/exe/Makenbi

$(LIBFEQIS): $(AWD)/src/feqis/*
	make -f $(AWD)/exe/Makefeqis lib

$(LIBSPIDER): $(AWD)/src/spider/*
	make -f $(AWD)/exe/Makespider
