import os, sys, logging
import const_text, eqns, eqns_init, config
import parse_as as pa

logger = logging.getLogger('as_parse.code_gen')
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

# Input:
#    tmp/model.tmp, astra_variables.json, ls fml/, ls fnc/
#
# Output:
#    src/tmp/*.f90, declar.fml, declar.fnc

awd = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


class CODE_GEN:


    def __init__(self, parse):


        eqns_lin = parse.eqns_lines

        linapp = False

#------------------------------- 
# Identifying SuBRoutines, SuBPprocesses
#------------------------------- 

        setv_sbr = ''
        detv_sbr = ''
        post_sbr = ''
        j_sbr = 1
        sbrs_d = {}

        for line in parse.sbr_lines:

            sbr_d = pa.parse_sbr(line)
            sbr_d['neq'] = j_sbr
            sbrs_d[line] = sbr_d
            locsbr = sbr_d['locsbr']
            a_str = ''
            if locsbr == -2:
                detv_sbr += pa.sbr_header(j_sbr, astr=a_str)
                detv_sbr += pa.write_sbr(sbr_d)
            elif locsbr == -1:
                detv_sbr += 'if (IPART == 2) then\n'
                detv_sbr += pa.sbr_header(j_sbr, astr=a_str)
                detv_sbr += pa.write_sbr(sbr_d)
                detv_sbr += 'endif\n'
            elif locsbr == 0:
                eqns_lin.append(line)
            elif locsbr == 1:
                post_sbr += pa.sbr_header(j_sbr, astr='')
                post_sbr += pa.write_sbr(sbr_d)

            j_arg = 1
            for key in ('dt', 'tmin', 'tmax', 'key'):
                val = sbr_d[key]
                if val != '':
                    try:
                        fval = float(val)
                        sval = pa.format_number(fval)
                        setv_sbr += 'DTEQ(%d, %d) = %s\n' %(j_arg, j_sbr, sval)
                    except:
                        l2f = pa.LINE2FOR(val, parse)
                        detv_sbr += 'DTEQ(%d, %d) = %s\n' %(j_arg, j_sbr, l2f)
                j_arg += 1
            j_sbr += 1

#-----------
# declar.fnc

        self.fnc = const_text.FNC.header + pa.write_declar(parse.fnc_list, ext='R')

#-----------
# declar.fml

        self.fml = pa.write_declar_fml(parse.fml_list)

#-----------
# associate_pointers

        self.associate_pointers = \
'''subroutine associate_pointers

use parameter_inc, only: NRD, n_sbr_max
use const_inc
use status_inc
use json_vars, only: n_const, n_var, n_varx, n_intern, n_intern2, n_prof, n_profx

implicit none

allocate(constValues(n_const))
allocate(varValues(n_var))
allocate(varxValues(n_varx))
allocate(internValues(n_intern + 4*n_sbr_max))
allocate(intern2Values(n_intern2))
allocate(profiles(NRD, n_prof))
allocate(profiles_x(NRD, n_profx))
allocate(DTEQ(4, n_sbr_max))

'''
        for j, const in enumerate(parse.constants):
            self.associate_pointers += '%s => constValues(%d)\n' %(const, j+1)
        self.associate_pointers += '\n'
        for j, var in enumerate(parse.variables):
            self.associate_pointers += '%s => varValues(%d)\n' %(var, j+1)
        self.associate_pointers += '\n'
        for j, varx in enumerate(parse.varx):
            self.associate_pointers += '%s => varxValues(%d)\n' %(varx, j+1)
        self.associate_pointers += '\n'
        for j, inter in enumerate(parse.intern1):
            self.associate_pointers += '%s => internValues(%d)\n' %(inter, j+1)
        self.associate_pointers += 'DTEQ(1:4, 1:n_sbr_max) => internValues(n_intern+1: n_intern + 4*n_sbr_max)\n\n'
        for j, inter in enumerate(parse.intern2):
            self.associate_pointers += '%s => intern2Values(%d)\n' %(inter, j+1)
        self.associate_pointers += '\n'
        for j, profx in enumerate(parse.profx):
            self.associate_pointers += '%s => profiles_x(:, %d)\n' %(profx, j+1)
        self.associate_pointers += '\n'
        for j, prof in enumerate(parse.prof):
            self.associate_pointers += '%s => profiles(:, %d)\n' %(prof, j+1)

        self.associate_pointers += \
'''
return
end subroutine associate_pointers
'''

#-----------
# postep.f90

        self.postep  = const_text.POSTEP.header
        self.postep += post_sbr
        self.postep += \
'''call markloc("tmp/postep.inc")

return
end subroutine POSTEP'''

#----------------------------
# detvar.f90

        detv_time = ''
        detv_rad  = ''

        for line in parse.detv_lines:
            line = line.strip()
            lbl = line.split('=', 1)[0].upper().strip()
            var = lbl.split('|', 1)[0] # For multiple statements with same left-hand side
            for jv, varm in enumerate(parse.variables):
                if var == varm:
                    jvar = jv + 1
                    detv_time += 'IFDFVX(%d) = max(IFDFVX(%d), 2)\n' %(jvar, jvar)
                    if '|' in line:
                        line = line.replace(lbl, var)
                    l2f = pa.LINE2FOR(line, parse)
                    detv_time += 'if (IFDFVX(%d) <= 2) %s\n'%(jvar, l2f)
                    break
            if var in parse.constants + parse.intern1 + parse.intern2:
                detv_time += pa.apptmp(lbl, parse)
            elif var in parse.profiles:
                detv_rad += pa.apptmp(lbl, parse)
            elif line[0] == '"' and line[-1] == '"':
                detv_rad += line[1: -1] + '\n'

        self.detvar  = const_text.DETVAR.header
        self.detvar += detv_time
        self.detvar += detv_sbr
        self.detvar += '''
! **** Radial profile computation
call markloc("detvar (profiles)")
do jdetv = 1, NA1
J = jdetv
'''
        self.detvar += detv_rad
        self.detvar += 'enddo\n'
        self.detvar += const_text.DETVAR.rad_tail
        self.detvar += \
'''
return
end subroutine DETVAR'''

#----------
# model.txt

        self.mtxt  = '''
 =====   Variables definition   =====
 =====   Radial profiles output   =====
  #  Scale  Name  Output expression
'''
        for jvar, var in enumerate(parse.asnamer):
            if var:
                out = var + '(r)'
            else:
                out = '0.d0'

            self.mtxt += '%3d  %6s %6s%s\n' %(jvar, parse.scaler[jvar].ljust(6), parse.namer[jvar].ljust(6), out)
        self.mtxt += ' =====   Time dependent values output   =====\n'
        self.mtxt += '  #  Scale  Name  Output expression\n'
        for jvar, var in enumerate(parse.asnamet):
            if var:
                tmp1 = var[:-1]
                out = var
                if var[-1] == 'B':
                    if tmp1 in parse.fnc_list + parse.profiles:
                        out = var.replace('B', '(a)')
                self.mtxt += '%3d  %6s %6s%s\n' %(jvar, parse.scalet[jvar].ljust(6), parse.namet[jvar].ljust(6), out)

#-----------
# inivar

        inivar = ''

# Dummy transport equations

        for jf in range(10):
            fj   = 'F%d' %jf
            inivar += 'do J=1, NA1\n'
            if fj in parse.var_defined:
                inivar += pa.apptmp(fj, parse)
            else:
                inivar += '%s(J) = 1.\n' %fj
            inivar += 'enddo\n'

        for lbl in ['NE', 'TE', 'TI', 'UPAR']:
# Initial condition
            if (parse.assign_d[lbl] == 'Missing') or (parse.assign_d[lbl] == 'AS' and parse.init_d[lbl] == ''):
                if lbl in ('TE', 'TI'):
                    inivar += pa.undef_inivar(lbl, config.short_d[lbl], defl='10 eV')
                if lbl == 'NE':
                    inivar += const_text.INIVAR.ne
                if lbl == 'UPAR':
                    inivar += pa.undef_inivar(lbl, config.short_d[lbl], defl='0 m/s', varx2='VTORX', varx3='VTORX(J)*IPOL(J)*RTOR/(RTOR+SHIF(J))')
                if parse.assign_d[lbl] != 'Missing':
                    logger.warning( 'Initial condition for %s is not defined' %lbl)
                    warn = '  %s=%sX(TSTART) will be used' %(lbl,lbl)
                    if lbl == 'UPAR':
                        warn.replace('UPARX', 'VTORX')
                    logger.warning(warn)
            else: # 'EQ', or 'AS' with explicit initial condition
                inivar += 'do J=1, NA1\n'
                inivar += pa.apptmp(lbl, parse)
                inivar += 'enddo\n'
            if parse.assign_d[lbl][:2] == 'EQ' and linapp:
                varb = config.bnd_d[lbl][0]
                short = config.short_d[lbl]
                rlbl = 'RO%s' %short 
                if varb in parse.init_d.keys():
                    inivar += 'ND1 = NA1\n'
                    inivar += pa.apptmp(varb, parse)
                if rlbl in parse.var_defined:
                    inivar += pa.apptmp(rlbl, parse)
                else:
                    inivar += 'ND1 = NA1\n'

        inivar += const_text.INIVAR.ni

# Current

        inivar += 'do J=1, NA1\n'
        if 'CC' in parse.var_defined:
            inivar += pa.apptmp('CC', parse)
        else:
            inivar += pa.LINE2FOR('CC = CCSP', parse)
        if 'CU' not in parse.var_defined:
            inivar += pa.apptmp('MU', parse)
            if parse.assign_d['CU'] == 'AS':
                inivar += 'CU(J) = CC(J)\n'
        inivar += 'enddo\n'

        self.iniv = inivar

        self.inivar  = const_text.INIVAR.header
        self.inivar += inivar
        self.inivar += \
'''
return
end subroutine INIVAR'''

#-----------
# ininam.f90

        inam  = 'AWD = "%s"\n' %awd
        if parse.arxuse:
            for j, arx in enumerate(parse.arxuse):
                inam += 'ARXUSE(%d) = %d\n' %(j+1, arx)
        for jlbl, lbl in enumerate(config.eqn_list):
            inam += 'LEQ(%d) = %d\n' %(jlbl+1, parse.leq_d[lbl])
        inam += const_text.ININAM.sb
        inam += 'call markloc("ininam")\n'
        inam += 'n_sbr = %d\n' %len(parse.sbr_lines)
        inam += 'NTOUT = %d\n' %len(parse.namet)
        inam += 'NROUT = %d\n' %len(parse.namer)
        inam += 'NXOUT = %d\n' %len(parse.namex)

        for jr, name in enumerate(parse.namer):
            if parse.scaler[jr] != '':
                inam += 'SCALER(%d) = %s\n' %(jr+1, parse.scaler[jr])
            inam += 'NAMER (%d) = "%s"\n' %(jr+1, name.ljust(4))
        for jrx, name in enumerate(parse.namex):
            inam += 'NAMEX (%d) = "%s"\n' %(jrx+1, name.ljust(6))
            inam += 'NWINDX(%d) = %s\n' %(jrx+1, parse.nwindx[jrx])
        for jt, name in enumerate(parse.namet):
            if parse.scalet[jt] != '':
                inam += 'SCALET(%d) = %s\n' %(jt+1, parse.scalet[jt])
            inam += 'NAMET (%d) = "%s"\n' %(jt+1, name.ljust(4))
        for line in parse.sbr_lines:
            j_sbr  = sbrs_d[line]['neq']
            sbrnam = sbrs_d[line]['name']
            inam += 'sbr_name(%d) = "%s"\n' %(j_sbr, sbrnam)
            inam += 'DTNAME(%d*4+n_intern) = "%s"//char(0)\n' %(j_sbr, sbrnam[:6])

        self.ininam  = const_text.ININAM.header
        self.ininam += inam
        self.ininam += \
'''
return
end subroutine ininam'''

#-----------
# setvar.f90

        setv_rho = ''
        for lbl in config.eqn_list:
            if lbl not in ('CU', 'Equil'):
                if parse.leq_d[lbl] == -1:
                    setv_rho += 'RO%s = ROC\n' %config.short_d[lbl]

        self.setvar = const_text.SETVAR.header
        self.setvar += setv_rho
        self.setvar += setv_sbr
        self.setvar += \
'''
return
end subroutine setvar'''

#--------------
# astra_out.f90
 
        self.astra_out  = const_text.RADOUT.header
        for jsgr, name in enumerate(parse.asnamer):
            if (jsgr > 128):
                raise ValueError('\n\n No more than 128 radial profiles in the output part (parameter NRW in astra_modules.f90)! \n\n') 
            if name:
                self.astra_out += pa.LINE2FOR('ROUT(J, %d) = %s' %(jsgr+1, name), parse)
            else:
                self.astra_out += 'ROUT(J, %d) = 0.d0\n' %(jsgr+1)
        self.astra_out += const_text.TIMOUT.header
        for jsig, name in enumerate(parse.asnamet):
            if (jsig > 128):
                raise ValueError('\n\n No more than 128 time traces in the output part (parameter NRW in astra_modules.f90)! \n\n') 
            if name:
                self.astra_out +=  pa.LINE2FOR('TOUT(LTOUT, %d) = %s' %(jsig+1, name), parse)
            else:
                self.astra_out += 'TOUT(LTOUT, %d) = 0.d0\n' %(jsig+1)
        self.astra_out += \
'''
return
end subroutine TIMOUT'''

#--------------------------------
# eqns_inc.f90, init_converge_step.f90

        eqns_tmp, init_tmp = eqns_init.eqns_init(parse)

        self.init_converge_step  = const_text.INIT_CONVERGE_STEP.header
        self.init_converge_step += init_tmp
        self.init_converge_step += const_text.INIT_CONVERGE_STEP.tail

        self.eqns_inc  = const_text.EQNS_INC.header
        self.eqns_inc += eqns_tmp
        self.eqns_inc += const_text.EQNS_INC.tail
