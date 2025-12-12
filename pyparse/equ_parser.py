import os, sys, logging, json
import const_text, eqns
import config
from parse_as import equ_prepare

logger = logging.getLogger('as_parse.equ_parser')
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)


class EQU_PARSER:


    def __init__(self, f_equ=None):

        self.f_equ = f_equ

#---------
# FML, FNC
#---------

        self.fml_list = []
        self.fnc_list = []

        for formula in sorted(os.listdir(config.fml_dir)):
            if os.path.isfile('%s/%s' %(config.fml_dir, formula)) and ('~' not in formula) and ('.' not in formula):
                self.fml_list.append(formula.upper())
        for fnc in sorted(os.listdir(config.fnc_dir)):
            if os.path.isfile('%s/%s' %(config.fnc_dir, fnc)) and (fnc[-1] in ('f','0','F','c', 'C')):
                tmp = fnc.split('.')[0]
                fnc = tmp.upper()
                self.fnc_list.append(fnc)

        if len(self.fml_list) > config.NFML:
            logger.error('Error, total amount of formulas > %d', config.NFML)
        if len(self.fnc_list) > config.NFML:
            logger.error('Error, total amount of functions > %d', config.NFML)

        for formula in self.fnc_list:
            if formula+'R' in self.fml_list:
                logger.error('Error! %s is both fml and fnc', formula)
                sys.exit()

#--------------------------------------
# Lists of arrays, constants, variables
#--------------------------------------

        f_json = '%s/astra_variables.json' %config.awd
        with open(f_json, 'r') as fjson:
            json_d = json.load(fjson)
        json_keys = {key: list(val.keys()) for key, val in json_d.items()}
        self.prof      = json_keys['profiles']
        self.profx     = json_keys['profiles_x']
        prof_ext = json_keys['strahl']
        self.profiles  = self.prof + self.profx + prof_ext
        self.constants = json_keys['constants']
        self.intern1   = json_keys['internal']
        self.intern2   = json_keys['intern2']
        self.variables = json_keys['variables']
        self.varx      = json_keys['variables_x']

#--------------------        
        self.arr_nam2 = []
        for arlist in config.coeff_d.values():
            self.arr_nam2 += arlist
        for arlist in config.flux_d.values():
            self.arr_nam2 += arlist

# Define short, initialise leq_d, assign_d

        self.leq_d = {}
        self.assign_d  = {}

        for eqn in config.eqn_list:
            self.leq_d[eqn] = -1 # No equation
            self.assign_d[eqn] = 'Missing'

#----------------------#
# Pre-parsing equ file #
#----------------------#
        if f_equ is None:
            return
        fequ_lin = equ_prepare(f_equ)

        equ_lines = []
        self.radout = [] # Lines with radial profiles output
        self.timout = [] # Lines with time traces output

        for line in fequ_lin:
            if ('=' in line) or (':' in line):
                equ_lines.append(line.replace(' ',''))
            if '\\' in line:
                self.radout.append(line)
            if '_' in line:
                flag_tim = True
                for sym in (':', '\\', '='):
                    if sym in line:
                        flag_tim = False
                        break
                if flag_tim:
                    self.timout.append(line)

#------------------------------
# Main parsing line by line, to preserve the sequence
#------------------------------

        self.detv_lines = []
        self.eqns_lines = []
        self.statements = []

        self.right_hand_d = {}
        self.right_hand_count = {}
        self.sbr_lines = []
        fluxes = []
        coeffs = []
        for flux in config.flux_d.values():
            fluxes += flux
        for coeff in config.coeff_d.values():
            coeffs += coeff

        for line in equ_lines:
            tmp = line.replace(':=', '=').replace('=>', '=').upper()

            if '=' in tmp:
                key, val = tmp.split('=', 1)
                if key in self.right_hand_d.keys(): # Variable was already defined above
                    self.right_hand_count[key] += 1
                    lbl = '%s|%d' %(key, self.right_hand_count[key]) # Keep several equ-lines with same left-hand side
                    self.right_hand_d[lbl] = val
                else:
                    self.right_hand_count[key] = 0
                    lbl = key
                self.right_hand_d[lbl] = val
                self.statements.append(line)
                if key in config.eqn_list + fluxes + coeffs:
                    self.eqns_lines.append(line)
                else:
                    self.detv_lines.append('%s = %s' %(lbl, val))

            if ':' in tmp:
# Equation
                key, val = tmp.split(':', 1)
                if key in config.eqn_list:
                    self.assign_d[key] = val
                    self.eqns_lines.append(line)
                else:
                    if key.endswith('*') and key[:-1] in ('TE', 'TI'):
                        self.assign_d[key[:-1]] = 'implicit_%s' %val
                        self.eqns_lines.append(line)
                    else:
                        self.sbr_lines.append(line)
        self.var_defined = [x for x in self.right_hand_d.keys() if '|' not in x]

# Equations, assignments

# LEQ
        #   -1  no equation
        #   0   type: AS
        #   1   type: EQ (default)
        #   2   type: FU
        #   3   type: heat conductivity flux + heat convection

        equ_type = None
        equ_solver = { \
            'EQ0':-2, # Ecy - Cylindrical
            'EQX':-1, # Ext - External file
            'EQD': 0, # Edf - Guessed equilibrium
            'E3M': 1, # EmEq - 3 moments
            'ESC': 2, # ESC
            'ESP': 3, # Spider
            'SPIDER': 3, # Spider
            'EQS': 4}  # SCoPE

        for key, val in self.assign_d.items():
            if val[:2] in ('', 'EQ'):
                self.leq_d[key] = 1 # EQ, default, such as TE:;
            if val == 'AS':
                self.leq_d[key] = 0
            elif val == 'FU':
                self.leq_d[key] = 2
        if equ_type is None:
            self.leq_d['Equil'] = -3
        else:
            self.leq_d['Equil'] = equ_solver[equ_type]

# Summary of assignments

        self.init_d = {}
        for var in config.eqn_list:
            self.init_d[var] = ''
            if self.assign_d[var] == 'Missing':
                self.init_d[var] = None
            else:
                for line in equ_lines:
                    linup = line.upper()
                    if var.upper() in linup and '=' in line:
                        tmp = line.split('=')
                        if var.upper().strip() == tmp[0].strip():
                            self.init_d[var] = tmp[1]
                logger.debug('%s %s %s %s', var, self.assign_d[var], self.init_d[var])

        arname = []
        equtxt = "".join(equ_lines)

        for j_eqn, eqn in enumerate(config.eqn_list):
            varx = eqn + 'X'
            if varx in equtxt:
                arname.append(varx)

        for lin in equ_lines:
            for jarr, var in enumerate(self.profx):
                if var in lin:
                    if var not in arname:
                        arname.append(var)
                        break

# Ordinal number in array list
        self.arxuse = []
        for varx in arname:
            for jarr, arr in enumerate(self.profx):
                if arr == varx:
                    self.arxuse.append(jarr+1)
                    break

# Time output

        self.namet   = []
        self.scalet  = []
        self.asnamet = []

        for sig in self.timout:
            tmp = sig.split('_')
            self.namet.append(tmp[0])     # used in ininam.tmp, model.txt
            self.asnamet.append(tmp[1].strip().upper())  # used in timout.tmp
            scale = ''
            if len(tmp) > 2:
                scale = tmp[2]
            self.scalet.append(scale)

# Profile output

        self.namer   = []
        self.asnamer = []
        self.scaler  = []
        self.namex   = []
        self.nwindx  = []

        for n_sgr, sgr in enumerate(self.radout):
            if n_sgr == config.NRW:
                logger.warning('Too many profiles for graphic output, cutting at 128')
                break
            tmp = sgr.split('\\')
            self.namer.append(tmp[0])
            self.asnamer.append(tmp[1].upper()) # not stored in ininam.tmp!
            scale = ''
            if '\\\\' in sgr:
                tmp = sgr.split('\\\\')[-1]
                if '\\' in tmp:
                    self.namex.append(tmp.split('\\')[0].upper())
                    scale = tmp.split('\\')[1]
                else:
                    self.namex.append(tmp)
                self.nwindx.append(n_sgr + 1)
            else:
                if len(tmp) > 2:
                    scale = tmp[2]
            self.scaler.append(scale)
