LIB_DIR = $(AWD)/lib
TSK_DIR = $(AWD)/bin
FOR_DIR = $(AWD)/main
MOD_DIR = $(AWD)/mod/for
TMP_DIR = $(AWD)/tmp
ATASK = $(TSK_DIR)/$(MODEL).exe
EQU  = $(AWD)/equ/$(MODEL)
LIBSBR = $(LIB_DIR)/libsbr.a
LIBFNC = $(LIB_DIR)/libfnc.a
LIBIPC = $(LIB_DIR)/libipc.a
LIBFOR = $(LIB_DIR)/libfor.a
LIBNBI = $(LIB_DIR)/libnbi.a
LIBFEQIS  = $(LIB_DIR)/libfeqis.a 
LIBSPIDER = $(LIB_DIR)/libspider.a 

MKL_LIB = -Wl,--start-group $(MKL_LIBDIR)/libmkl_intel_lp64.a $(MKL_LIBDIR)/libmkl_core.a $(MKL_LIBDIR)/libmkl_sequential.a -Wl,--end-group
LIB = -lnbi -lipc -lsbr -lfnc $(XLBR) $(NCDF_LIB) $(QLKNN_LIB) $(TGLF_LIB) $(TORB_LIB) $(RABBIT_LIB) -lfor -lfeqis -lspider -L$(LIB_DIR) $(MKL_LIB)

MODTMP = $(AWD)/tmp/model.tmp
AENV = $(AWD)/exe/astra_rc
FOR_SRC = \
$(FOR_DIR)/astra_main.f90 \
$(FOR_DIR)/stepup.f90
TMP_SRC = \
$(TMP_DIR)/astra_out.f90 \
$(TMP_DIR)/inivar.f90 \
$(TMP_DIR)/setvar.f90 \
$(TMP_DIR)/detvar_init.f90 \
$(TMP_DIR)/detvar.f90 \
$(TMP_DIR)/postep.f90 \
$(TMP_DIR)/subproc.f90 \
$(TMP_DIR)/ininam.f90 \
$(TMP_DIR)/converge_init.f90 \
$(TMP_DIR)/eqns_inc.f90

PY_TMP_GEN = $(AWD)/pyparse/parser_main.py -awd $(AWD)

FMLS    := $(wildcard $(AWD)/fml/*)
PYFILES := $(wildcard $(AWD)/pyparse/*.py)

# Make Commands
.PHONY: all
all: directories update_model $(ATASK)

.PHONY: directories
directories:
	@mkdir -p $(TSK_DIR)
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(AWD)/.res/ncdf
	@mkdir -p $(AWD)/dat

.PHONY: update_model
update_model:
	rsync -c $(EQU) $(MODTMP)

# Linking, create ASTRA executable

$(ATASK): $(LIBFOR) $(LIBNBI) $(LIBFEQIS) $(LIBSPIDER) $(LIBSBR) $(LIBFNC) $(LIBIPC) $(AENV) $(MODTMP) $(FOR_SRC) $(PYFILES)
	$(PY_TMP_GEN)
	$(AFC) -qopenmp -g -traceback -o $(ATASK) $(FOR_SRC) $(TMP_SRC) -I$(MOD_DIR) -I$(NBI_INC) $(LIB)

# Generate FORTRAN code
# EQU is copied into MODTMP above if they are different, here if EQU is touched
$(MODTMP): $(EQU)
	cp -p $(EQU) $(MODTMP)

# User defined sbr and fnc

$(LIBSBR): $(AWD)/sbr/*
	make -f $(AWD)/exe/Makesub sub=sbr

$(LIBFNC): $(AWD)/fnc/*.f90
	make -f $(AWD)/exe/Makesub sub=fnc

$(LIBIPC): $(AWD)/ipc/*
	make -f $(AWD)/exe/Makesub sub=ipc

# Former 'server' part

$(LIBFOR): $(AWD)/src/for/*
	make -f $(AWD)/exe/Makefor

$(LIBNBI): $(AWD)/src/nbi/*
	cd $(AWD)/src; make -f $(AWD)/exe/Makenbi

$(LIBFEQIS): $(AWD)/src/feqis/*
	make -f $(AWD)/exe/Makefeqis

$(LIBSPIDER): $(AWD)/src/spider/*
	make -f $(AWD)/exe/Makespider
