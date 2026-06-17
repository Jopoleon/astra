import logging
import const_text, config
import parse_as as pa

logger = logging.getLogger('as_parse.eqns')
logger.setLevel(logging.INFO)

# var_defined: list of all variables assigned in the equ file (left hand side of '=' statements)
# assign_type: 'Missing', 'AS', 'EQ', 'EQ[j,rho]', 'FU'?

def none_in(list_in, var_list):
    return not count_in(list_in, var_list)

def count_in(list_in, var_list):
    return sum((lbl in var_list) for lbl in list_in)

def pre_eqn(parse, key, assign_type=None):

    if assign_type is None:
        assign_type = parse.assign_d[key]

    var_defined = [x for x in parse.var_defined]

    pre_txt = ''
    if assign_type in ('AS', 'Missing'):
        pre_txt += '! **** %s assignment\n' %config.labels_d[key]
        pre_txt += 'call markloc("%s assignment")\n' %key
    elif assign_type[:2] == 'EQ':
        if none_in(config.coeff_d[key] + config.flux_d[key], var_defined):
            for varb in config.bnd_d[key]:
                if varb in var_defined:
                    logger.warning('Equation for "%s" requested, but not defined', key)
                    logger.warning('Boundary condition %s = %s ignored', varb, parse.right_hand_d[varb])
                    var_defined.remove(varb)
        pre_txt += '! **** %s equation\n' %config.labels_d[key]
        pre_txt += 'call markloc("%s equation")\n' %key

    pre_txt += 'do J=1, NA1\n'
    for var in config.coeff_d[key]:
        pre_txt += pa.apptmp(var, parse)
    for var in config.flux_d[key]:
        if var in var_defined:
            pre_txt += pa.apptmp(var, parse)
        else:
            pre_txt += '%s(J) = 0.\n' %var
    if key != 'UPAR':
        flux = config.flux_d[key][0]
        pre_txt += '%sTOT(J) = %s(J)\n' %(flux, flux)
        pre_txt += 'enddo\n'

    return pre_txt, var_defined


def rhoBC(assign_type):

    rho_val = None
    tbeg = -1.e12
    tend = 1.e12
    if '[' in assign_type:
        tmp = assign_type.split('[')[1].split(']')[0]
        pieces = tmp.split(',')
        n_commas = len(pieces) - 1
        if n_commas == 0:
            rho_val = 'RFA(%s)' %pieces[0]
        else:
            jprof = int(tmp[0])
            if jprof == 1:
                rho_val = 'RFAN(%s)' %pieces[1]
            elif jprof == 2:
                rho_val = '%s*ROC' %pieces[1]
            else:
                logger.error('First argument in ...:EQ[ , ] can be only either 1 or 2')
                logger.error('Please amend your equ file')
                sys.exit()
            if n_commas == 2:
                tbeg = float(pieces[2])
            elif n_commas == 3:
                tbeg = float(pieces[2])
                tend = float(pieces[3])

    return rho_val, tbeg, tend


def linInterp2sep(var):
    return '''
if (ND1 < NA) then
YA = (%s(NA1) - %s(ND1))/(ROC - RHO(ND1))
YB = (%s(ND1)*ROC - %s(NA1)*RHO(ND1))/(ROC - RHO(ND1))
do j=ND1+1, NA
%s(j) = YB + RHO(j)*YA
enddo
endif
'''  %(var, var, var, var, var)


def bnd_init(var, var_defined, parse):
    '''Writing boundary/external condition for TE, TI, NE, F*, UPAR equations'''

    assign_type = parse.assign_d[var]
    bnd_txt = ''
    rho_bnd, tbeg, tend = rhoBC(assign_type)
    if rho_bnd is None:
        bnd_txt += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        bnd_txt += 'tbeg_eq = %12.4e\n' %tbeg
        bnd_txt += 'tend_eq = %12.4e\n' %tend
        bnd_txt += 'if (TIME >= tbeg_eq .and. TIME <= tend_eq) then\n'
        bnd_txt += 'ND1 = NODE(%s)\n' %l2f.strip()
        if var in var_defined:
            bnd_txt += 'do j=ND1, NA1\n'
            bnd_txt += pa.apptmp(var, parse)
            bnd_txt += 'enddo\n'
        bnd_txt += 'else if (TIME < tbeg_eq) then\n'
        bnd_txt += 'ND1 = 1\n'
        if var in var_defined:
            bnd_txt += 'do j=ND1, NA1\n'
            bnd_txt += pa.apptmp(var, parse)
            bnd_txt += 'enddo\n'
        bnd_txt += 'else if (TIME > tend_eq) then\n'
        bnd_txt += 'ND1 = 1\n'
        if var in var_defined:
            bnd_txt += 'do j=ND1, NA1\n'
            bnd_txt += '%s(j) = %sO(j)\n' %(var, var)
            bnd_txt += 'enddo\n'
        bnd_txt += 'endif\n'
    if var == 'UPAR':
        varx = 'VTORX'
    else:
        varx = '%sX' %var
    varb_list = config.bnd_d[var]
    bnd_count = count_in(varb_list, var_defined)
    if bnd_count == 0: # No boundary condition is set
        if var not in var_defined:
            logger.warning('Neither BC nor external condition for %s is given', var)
            logger.warning('Assuming %s(ND1: NA1)=%s(ND1: NA1)', var, varx)
            bnd_txt += '%s(ND1: NA1) = %s(ND1: NA1)\n' %(var, varx)
        bnd_txt += '%sO(ND1: NA1) = %s(ND1: NA1)\n' %(var, var)
        bnd_txt += 'bctype = 1\n'
    elif bnd_count == 1:
        if varb_list[0] in var_defined: # NEB
            bnd_txt += '%s(ND1) = %s' %(var, pa.LINE2FOR(parse.right_hand_d[varb_list[0]], parse) )
            bnd_txt += '%sO(ND1: NA1) = %s(ND1: NA1)\n' %(var, var)
            bnd_txt += 'bctype = 1\n'
            if var not in var_defined: # Linear extrapolation towards rho=1
                logger.warning('BC defined, but no external condition for %s', var)
                logger.warning('Assuming linear interpolation between rhoBC and separatrix')
                bnd_txt += linInterp2sep(var)
        else: # QNB
            flux = varb_list[1][:-1]     # flux is 'QN', 'QE',...
            if varb_list[1] in var_defined:
                bnd_txt += '%s(ND1) = %s' %(flux, pa.LINE2FOR(parse.right_hand_d[varb_list[1]], parse) )
            else:
                bnd_txt += '%s(ND1) = %s(ND1)*(%s)\n' %(flux, var, pa.LINE2FOR(parse.right_hand_d[varb_list[2]], parse).strip() )
            bnd_txt += 'bc_values(1) = %s(ND1)\n' %flux
            bnd_txt += 'bctype = 2\n'
    else:
        raise ValueError('Too many boundary conditions defined for %s in the equ file', var)

    return bnd_txt


def cuasn(parse, bc='CU', neq=1):

    cuas_txt = const_text.CUAS.header
 
    if 'MV' in parse.var_defined:
         logger.warning('MV is not used to define FV')
    else:
        cuas_txt += 'FV(1: NA1) = 0.\n'

    cuas_txt += 'do j=1, NA1\n'
    for var in ('DC', 'HC', 'XC', 'CD', 'CC'):
        cuas_txt += pa.apptmp(var, parse)

    cuas_txt += cubs(parse.var_defined)

    if 'CD' in parse.var_defined:
        cuas_txt += 'YWA(J) = CUBS(J) + CD(J)\n'
    else:
        cuas_txt += 'YWA(J) = CUBS(J)\n'

    cuas_txt += 'enddo ! j (radial loop)\n'
    if neq == 0:
        cuas_txt += 'FP(1) = FV(1)\n'
    else:
        cuas_txt += 'FP(1) = FPO(1) + UPL(1)*TAU\n'

    if bc == 'MU': # MU
        cuas_txt += 'do j=1, NA1\n'
        cuas_txt += pa.apptmp('MU', parse)
        cuas_txt += const_text.CUAS.mu
    else:
        cuas_txt += const_text.CUAS.cu1
        cuas_txt += pa.apptmp('CU', parse)
        cuas_txt += const_text.CUAS.cu2
        if 'MV' in parse.var_defined:
            cuas_txt += 'MU(J) = YJ_CU*MU(J) + MV(j)\n'
            cuas_txt += 'FP(J) = YJ_CU*(FP(J) - FP(1)) + FP(1) + FV(J)\n'
        else:
            cuas_txt += 'MU(J) = YJ_CU*MU(J)\n'
            cuas_txt += 'FP(J) = YJ_CU*(FP(J) - FP(1)) + FP(1)\n'

    cuas_txt += const_text.CUAS.cu_mu

    return cuas_txt


def cubs(var_defined):

    cubs_line = ''
    if none_in(['DC', 'HC', 'XC'], var_defined):
        if 'CUBS' in var_defined:
            cubs_line += pa.apptmp('CUBS', parse)
        else:
            cubs_line += 'CUBS(J) = 0.\n'
    else:
        cubs_line += '''if (j == NA1) then
CUBS(NA1) = max(0., 2.*CUBS(NA) - CUBS(NA-1))
else
CUBS(J) = YA*(FP(J+1) - FP(J))*( 0.'''
        if 'HC' in var_defined:
            cubs_line += ' + HC(J)*(TE(J+1) - TE(J))/(TE(J+1) + TE(J))'
        if 'XC' in var_defined:
            cubs_line += ' + XC(J)*(TI(J+1) - TI(J))/(TI(J+1) + TI(J))'
        if 'DC' in var_defined:
            cubs_line += ' + DC(J)*(NE(J+1) - NE(J))/(NE(J+1) + NE(J))'
        cubs_line += ' ) * 0.5*(IPOL(J) + IPOL(J+1))\n'
        cubs_line += 'endif\n'

    return cubs_line


def neeqn(parse):

    assign_type = parse.assign_d['NE']

    ne_txt, var_defined = pre_eqn(parse, 'NE')
    if assign_type == 'Missing':
        return ne_txt

    ne_txt += 'do J=1, NA\n'
    ne_txt += 'YWA(J) = 0.'
    if 'DN' in var_defined:
        ne_txt += ' + DN(J)'
    if 'DVN' in var_defined:
        ne_txt += ' + DVN(J)'
    ne_txt += '\n'
    if 'DVN' in var_defined:
        ne_txt += 'YWB(J) = DVN(j)*log(NE(j)/NE(j+1))/HRO'
    else:
        ne_txt += 'YWB(J) = 0.'
    if 'HN' in var_defined:
        ne_txt += ' + 2.*HN(J)/(TE(J+1) + TE(J))*(TE(J+1) - TE(J))/HRO'
    if 'XN' in var_defined:
        ne_txt += ' + 2.*XN(J)/(TI(J+1) + TI(J))*(TI(J+1) - TI(J))/HRO'
    if 'CN' in var_defined:
        ne_txt += ' - CN(J)'
    ne_txt += '\n'
    ne_txt += 'enddo\n'

    if assign_type == 'AS':
        ne_txt += 'do j=1, NA1\n'
        if 'NE' in var_defined:
            ne_txt += pa.apptmp('NE', parse)
        else:
            logger.warning('Missing NE assignment, assuming NEX(t)')
            ne_txt += 'NE(J) = NEX(J)\n'
        ne_txt += 'enddo\n'
        ne_txt += const_text.NEEQN.assigned
    elif assign_type[:2] == 'EQ':
        ne_txt += bnd_init('NE', var_defined, parse)
        ne_txt += const_text.NEEQN.eqn

    if 'SNN' in var_defined:
        ne_txt += 'SNTOT(1: NA) = SNTOT(1: NA) + SNN(1: NA)*NE(1: NA)\n'
    ne_txt += 'GN(NA1) = QN(NA1)/SLAT(NA1)\n'
    ne_txt += 'SNTOT(NA1) = SNTOT(NA)\n\n'

    return ne_txt


def tieqn(parse):

    assign_type = parse.assign_d['TI']

    ti_txt, var_defined = pre_eqn(parse, 'TI')
    if assign_type == 'Missing':
        return ti_txt

    ti_txt += 'do J=1, NA\n'
    ti_txt += 'YWA(J) = 0.'
    if 'XI' in var_defined:
        ti_txt += ' + XI(J)' 
    if 'DVI' in var_defined:
        ti_txt += ' + DVI(J)'
    if 'XN' in var_defined and assign_type[:2] == 'EQ':
        ti_txt += ' + GN2I*XN(J)'
    ti_txt += '\n'
    ti_txt += 'YWA(J) = YWA(J)*(NI(J+1) + NI(J))*0.5\n'
    LB = True
    LC = True
    if 'DN' in var_defined and assign_type[:2] == 'EQ':
        if 'DI' in var_defined:
            ti_txt += 'YWB(J) = (DI(J) + GN2I*DN(J))*(NI(J+1) + NI(J))/(NE(J+1) + NE(J))\n'
        else:
            ti_txt += 'YWB(J) = GN2I*DN(J)*(NI(J+1) + NI(J))/(NE(J+1) + NE(J))\n'
    else:
        if 'DI' in var_defined:
            ti_txt += 'YWB(J) = DI(J)*(NI(J+1) + NI(J))/(NE(J+1) + NE(J))\n'
        else:
            LB = False
    if 'HN' in var_defined and assign_type[:2] == 'EQ':
        if 'HI' in var_defined:
            ti_txt += 'YWC(J) = (HI(J) + GN2I*HN(J))*(NI(J+1) + NI(J))/(TE(J+1) + TE(J))\n'
        else:
            ti_txt += 'YWC(J) = GN2I*HN(J)*(NI(J+1) + NI(J))/(TE(J+1) + TE(J))\n'
    else:
        if 'HI' in var_defined:
            ti_txt += 'YWC(J) = HI(J)*(NI(J+1) + NI(J))/(TE(J+1) + TE(J))\n'
        else:
            LC = False
    ti_txt += 'YWD(J) = 2.*GN2I*GNX(J)*SLAT(J)/G11(J)/(NE(J+1) + NE(J))'
    if 'CI' in var_defined:
        ti_txt += ' + CI(J)'
    if 'CN' in var_defined:
        ti_txt += ' + GN2I*CN(J)'
    ti_txt += '\n'
    ti_txt += 'YWD(J) = 0.5*YWD(J)*(NI(J+1) + NI(J))\n'
    ti_txt += 'YWB(J) = -YWD(J)'
    if LB:
        ti_txt += ' + YWB(J)*(NE(J+1) - NE(J))/HRO'
    if LC:
        ti_txt += ' + YWC(J)*(TE(J+1) - TE(J))/HRO'
    ti_txt += '\n'
    if 'DVI' in var_defined:
        ti_txt += 'YWB(J) = YWB(J) + DVI(J)*0.5*(NI(j) + NI(j+1))*log(TI(j)/TI(j+1))/HRO\n'
    ti_txt += 'enddo\n'

    if assign_type == 'AS':
        ti_txt += 'do j=1, NA1\n'
        if 'TI' in var_defined:
            ti_txt += pa.apptmp('TI', parse)
        else:
            logger.warning('Missing TI assignment, assuming TIX(t)')
            ti_txt += 'TI(J) = TIX(J)\n'
        ti_txt += 'enddo\n'
        ti_txt += const_text.TIEQN.assigned
    elif assign_type[:2] == 'EQ':
        ti_txt += bnd_init('TI', var_defined, parse)
        if 'DSI' in var_defined:
            ti_txt += 'DSI(ND1) = 1.\n'
        else:
            ti_txt += 'DSI(ND1) = 0.\n'
        ti_txt += const_text.TIEQN.eqn

    return ti_txt


def teeqn(parse):

    assign_type = parse.assign_d['TE']

    te_txt, var_defined = pre_eqn(parse, 'TE')
    if assign_type == 'Missing':
        return te_txt

    te_txt += 'do J=1, NA\n'
    te_txt += 'YWA(J) = 0.'
    if 'HE' in var_defined:
        te_txt += ' + HE(J)' 
    if 'DVE' in var_defined:
        te_txt += ' + DVE(J)'
    if 'HN' in var_defined and assign_type[:2] == 'EQ':
        te_txt += ' + GN2E*HN(J)'
    te_txt += '\n'
    te_txt += 'YWA(J) = YWA(J)*(NE(J+1) + NE(J))*0.5\n'
    LB = True
    LC = True
    if 'DN' in var_defined and assign_type[:2] == 'EQ':
        if 'DE' in var_defined:
            te_txt += 'YWB(J) = DE(J) + GN2E*DN(J)\n'
        else:
            te_txt += 'YWB(J) = GN2E*DN(J)\n'
    else:
        if 'DE' in var_defined:
            te_txt += 'YWB(J) = DE(J)\n'
        else:
            LB = False
    if 'XN' in var_defined and assign_type[:2] == 'EQ':
        if 'XE' in var_defined:
            te_txt += 'YWC(J) = (XE(J) + GN2E*XN(J))*(NE(J+1) + NE(J))/(TI(J+1) + TI(J))\n'
        else:
            te_txt += 'YWC(J) = GN2E*XN(J)*(NE(J+1) + NE(J))/(TI(J+1) + TI(J))\n'
    else:
        if 'XE' in var_defined:
            te_txt += 'YWC(J) = XE(J)*(NE(J+1) + NE(J))/(TI(J+1) + TI(J))\n'
        else:
            LC = False
    te_txt += 'YWD(J) = 0.'
    if 'CE' in var_defined:
        te_txt += ' + CE(J)'
    if 'CN' in var_defined:
        te_txt += ' + GN2E*CN(J)'
    te_txt += '\n'

    te_txt += 'YWD(J) = YWD(J)*(NE(J+1) + NE(J))*0.5 + GN2E*GNX(J)*SLAT(J)/G11(J)\n'
    te_txt += 'YWB(J) = -YWD(J)'
    if LB:
        te_txt += ' + YWB(J)*(NE(J+1) - NE(J))/HRO'
    if LC:
        te_txt += ' + YWC(J)*(TI(J+1) - TI(J))/HRO'
    te_txt += '\n'
    
    if 'DVE' in var_defined:
        te_txt +=  'YWB(J) = YWB(J) + DVE(J)*0.5*(NE(j) + NE(j+1))*log(TE(j)/TE(j+1))/HRO\n'
    te_txt += 'enddo\n'

    if assign_type == 'AS':
        te_txt += 'do j=1, NA1\n'
        if 'TE' in var_defined:
            te_txt += pa.apptmp('TE', parse)
        else:
            logger.warning('Missing TE assignment, assuming TEX(t)')
            te_txt += 'TE(J) = TEX(J)\n'
        te_txt += 'enddo\n'
        te_txt += const_text.TEEQN.assigned
    elif assign_type[:2] == 'EQ':
        te_txt += bnd_init('TE', var_defined, parse)
        if 'DSE' in var_defined:
            te_txt += 'DSE(ND1) = 1.\n'
        else:
            te_txt += 'DSE(ND1) = 0.\n'
        te_txt += const_text.TEEQN.eqn

    return te_txt


def cueqn(parse):

    cueq_txt = const_text.CUEQN.header
    if 'MV' in parse.var_defined:
         logger.warning('MV is not used to define FV')
         cueq_txt += pa.apptmp('MV', parse)
         cueq_txt += const_text.CUEQN.mv

    for coeff in ('DC', 'HC', 'XC', 'CD', 'CC'):
        cueq_txt += pa.apptmp(coeff, parse)

    cueq_txt += cubs(parse.var_defined)

    if 'CD' in parse.var_defined:
        cueq_txt += pa.apptmp('CD', parse)
        cueq_txt += 'YWD(J) = CUBS(J) + CD(J)\n'
    else:
        cueq_txt += 'YWD(J) = CUBS(J)\n'

    cueq_txt += '''
YWD(J) = YWD(J)*YD/(IPOL(J)**3 * G33(J))
YWB(J) = 0.4*GP*CC(J)*YC/IPOL(J)**2
enddo
'''

    if 'UEXT' not in parse.var_defined and 'LEXT' not in parse.var_defined and 'IPL' in parse.var_defined:
        cueq_txt += const_text.CUEQN.prescribed_ipl
    if 'UEXT' in parse.var_defined and 'LEXT' not in parse.var_defined:
        cueq_txt += const_text.CUEQN.prescribed_uloop
    if 'LEXT' in parse.var_defined:
        cueq_txt += const_text.CUEQN.circuit_eqn

    cueq_txt += const_text.CUEQN.eqn

    if 'LEXT' in parse.var_defined:
        cueq_txt += 'PSIFB = PSIEXT - PSPLEX*dfpdrbm12\n'

    cueq_txt += \
'''do J=1, NA1
UPL(J) = (FP(J) - FPO(J))/TAU - YQDCMF(J)
! this is dpsi/dt_rho. to get dpsi/dt_x full
!(so boundary values is dpsi/dt_x=1), one needs to add yqdcmf. Or simply compute it as dFP(NA1)/dt.
ULON(J) = IPOL(J)*G33(J)*(UPL(J) - GP2*ROC**2 * BTOR*BBDOT*MU(J)) !this is correct, also 
!goes into Ohmic power. Ohmic power is not computed with dPsi/dt_x, but dPsi/dt_phi (ULON)
enddo
call CUOFP()
'''

    if 'UEXT' in parse.var_defined or 'LEXT' in parse.var_defined:
        cueq_txt += '!DFPDR = ((FP(NA1) - FP(NA) - (FV(NA1) - FV(NA)))/HRO)\n'
        cueq_txt += 'IPL = 5.*IPOL(NA1)*G22(NA)*((FP(NA1) - FP(NA) - (FV(NA1) - FV(NA)))/HRO)/GP2/RTOR\n'

    return cueq_txt


def fjeqn(parse, jeq):

    key  = 'F%d'   %jeq
    ro   = 'RO%d'  %jeq
    df   = 'D%s'   %key
    dvf  = 'DV%s'  %key
    vf   = 'V%s'   %key
    qf   = 'Q%s'   %key
    qff  = 'QF%s'  %key
    sf   = 'S%s'   %key
    sff  = 'SF%s'  %key
    gf   = 'G%s'   %key
    varb = '%sB'   %key
    qfb  = 'Q%sB'  %key
    qffb = 'QF%sB' %key

    assign_type = parse.assign_d[key]
    if assign_type == 'Missing': # is actually prevented ahead
        return ''

    fj_txt, var_defined = pre_eqn(parse, key)

    fj_txt += 'do J=1, NA1\n'
    fj_txt += 'YWA(J) = 0.'
    if df in var_defined:
        fj_txt += ' + %s(J)' %df
    if dvf in var_defined:
        fj_txt += ' + %s(J)' %dvf
    fj_txt += '\n'
    fj_txt += 'enddo\n'

    fj_txt += 'do J=1, NA\n'
    if dvf in var_defined:
        lin2 = 'YWB(J) = %s(J)*log(%s(J)/%s(J+1))/HRO' %(dvf, key, key)
    else:
        lin2 = 'YWB(J) = 0.'
    if vf in var_defined:
        lin2 += '-%s(J)' %vf
    fj_txt += lin2 + '\n'
    fj_txt += 'enddo\n'

    if assign_type == 'AS':
        if key in var_defined:
            fj_txt += 'do j=1, NA1\n'
            fj_txt += pa.apptmp(key, parse)
            fj_txt += 'enddo\n'
        fj_txt += 'do j=1, NA\n'
        fj_txt += '%s(J) = -G11(J)*(YWA(J)*(%s(J+1) - %s(J))/HRO + 0.5*YWB(J)*(%s(J+1) + %s(J)))\n' %(qf, key, key, key, key)
        if gf in var_defined:
            fj_txt += '%s(J) = %s(J) + SLAT(J)*%s(J)\n' %(qf, qf, gf)
        fj_txt += 'enddo\n'
    else:
        fj_txt += bnd_init(key, var_defined, parse)
        fj_txt += const_text.FJEQN.eqn
        fj_txt += 'NA1%d = ND1\n' %jeq
        fj_txt += \
'call RUNEQ(YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), F%dO(1: NA1), YWWB(1: NA1), YVR(1: NA1), unit_coeff, G11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), SFF%d(1: NA1),SF%d(1: NA1), RBDOT, BBDOT, ND1, NA1, HRO, TAU, RHO(1: NA1), imethod, bctype, bc_values, F%d(1: NA1), QF%d(1: NA1), YQDCM(1: NA1), MPHIT(1: NA1))\n' %(jeq, jeq, jeq, jeq, jeq)

    if sff in var_defined:
        fj_txt += '%sTOT(NA1) = %sTOT(NA1) + %s(NA1)*%s(NA1)\n' %(sf, sf, sff, key)

    if qfb in var_defined and none_in([qffb, varb, ro], var_defined):
        fj_txt += '%s(NA1) = %sB\n' %(qf, qf)
    else:
        fj_txt += '%s(NA1) = %s(NA)\n' %(qf, qf)

    if gf not in var_defined:
        fj_txt += '%s(NA1) = %s(NA1)/SLAT(NA1)\n' %(gf, qf)

    fj_txt += '%sTOT(NA1) = %sTOT(NA)\n' %(sf, sf)

    return fj_txt


def upeqn(parse):

    assign_type = parse.assign_d['UPAR']
    if assign_type == 'Missing': # never occurs, prevented ahe
        return ''

    up_txt, var_defined = pre_eqn(parse, 'UPAR')

    if 'XUPAR' in var_defined:
        up_txt += 'YWA(J) = XUPAR(J)*MRHO(J)/IPOL(J)*RTOR\n'
    else:
        up_txt += 'YWA(J) = 0\n'

    up_txt += 'YWR(J) = 0.'
# Parallel stress tensor
# defined as residuai stress/grad rho**2
    if 'RUPAR' in var_defined:
        up_txt += ' + RUPAR(J)/G11(J)*VRS(J)'
    if 'RUPYR' in var_defined:
        up_txt += ' - RUPYR(J)/G11(J)*VRS(J)'
    if 'RUPFR' in var_defined:
        up_txt += ' + RUPFR(J)/G11(J)*VRS(J)'
    up_txt += '\n'
    up_txt += 'if (j <= NA) then\n'
    up_txt += 'YWD(J) = 0.\n'
    if 'CNPAR' in var_defined:
        up_txt += 'YWgradF(J) = 2.*(IPOL(J+1) - IPOL(J))/(IPOL(J) + IPOL(J+1))/HRO\n'
        up_txt += 'YWD(J) = YWD(J) + (CNPAR(J) + XUPAR(J)*YWgradF(J))*MRHO(J)/IPOL(J)*RTOR\n'
    if assign_type[:2] == 'EQ':
        up_txt += 'YWB(J) = -YWD(J)\n'   #pinch term
        up_txt += 'YWC(J) = 0.0*TTRQ(J)\n'

    up_txt += '''endif
UPS0(j) = MRHO(j)*G41(j)*RTOR/IPOL(j)
! Neoclassical corrections
UPS1(j) = -DLNEO(j)*G41(j)*RTOR/IPOL(j)
UPS2(j) = SGNEO(J)*RTOR*BTOR*IPOL(j)*(1. - BDB02(j)*G41(j)/IPOL(j)**0.2)
enddo
'''

    if assign_type == 'AS':
        up_txt += 'do j=1, NA1\n'
        if 'UPAR' in var_defined:
            up_txt += pa.apptmp('UPAR', parse)
        else:
            up_txt += 'UPAR(J) = VTORX(J)\n'
        up_txt += 'enddo\n'
        up_txt += const_text.UPEQN.assigned
    else:
        up_txt += bnd_init('UPAR', var_defined, parse)
        up_txt += const_text.UPEQN.eqn

    return up_txt


def tetieqn(parse):

    var_defined = parse.var_defined

    assign_type = parse.assign_d['TE']
    impl, asstyp = assign_type.split('_', 1)

    teti = ''
    if none_in(['DE', 'HE', 'XE', 'CE', 'PE', 'PET'], var_defined):
        for lbl in ('TEB', 'QEB', 'QETB'):
            if lbl in var_defined:
                logger.warning('Equation for "TE" requested but not defined')
                logger.warning('Boundary condition ignored')
                print(pa.apptmp(lbl, parse))
                var_defined.remove(lbl)

    if asstyp[:2] == 'EQ':
        teti += '! **** Electron temperature equation\n'
        teti += 'call markloc("TE equation")\n'
    else:
        logger.error('ERROR tetieqn: TE cannot be ASSIGNED with implicit scheme')
        return ''

    teti += 'do j=1, NA1\n'
    for lbl in ['DE', 'HE', 'XE', 'CE', 'PE', 'PET', 'DVE', 'DSE']:
        if lbl in var_defined:
            teti += pa.apptmp(lbl, parse)
    teti += 'if (j > NA) CYCLE\n'
    txt = 'YWA1(J) = 0.'
    if 'HE' in var_defined:
        txt += ' + HE(J)'
    if 'HN' in var_defined:
        txt += ' + GN2E*HN(J)'
    if 'DVE' in var_defined:
        txt += ' + DVE(J)'
    teti += txt + '\n'
    teti += 'YWA1(J) = YWA1(J)*(NE(J+1) + NE(J))*0.5\n'

    LB = 1
    if 'DN' in var_defined:
        if 'DE' in var_defined:
            teti += 'YWB1(J) = (DE(J) + GN2E*DN(J))\n'
        else:
            teti += 'YWB1(J) = GN2E*DN(J)\n'
    else:
        if 'DE' in var_defined:
            teti += 'YWB1(J) = DE(J)\n'
        else:
            LB = 0

    LC = 1
    if 'XN' in var_defined:
        if 'XE' in var_defined:
            teti += 'YWC(J) = (XE(J) + GN2E*XN(J))*(NE(J+1) + NE(J))/(TI(J+1) + TI(J))\n'
        else:
            teti += 'YWC(J) = GN2E*XN(J)*(NE(J+1) + NE(J))/(TI(J+1) + TI(J))\n'
    else:
        if 'XE' in var_defined:
            teti += 'YWC(J) = XE(J)*(NE(J+1) + NE(J))/(TI(J+1) + TI(J))\n'
        else:
            LC = 0

    txt = 'YWD(J) = 0.'
    if 'CE' in var_defined:
        txt += ' + CE(J)'
    if 'CN' in var_defined:
        txt += '+ GN2E*CN(J)'
    teti += txt + '\n'
    teti += 'YWD(J) = YWD(J)*(NE(J+1) + NE(J))*0.5 + GN2E*GNX(J)*SLAT(J)/G11(J)\n'
    teti += 'enddo\n'
    teti += 'do J=1, NA1\n'
    if 'PET' not in var_defined:
        teti += 'PET(J) = 0.\n'
    if 'PE' not in var_defined:
        teti += 'PE(J) = 0.\n'
    teti += 'PETOT(J) = PE(J)\n'
    teti += 'if (j > NA) CYCLE\n'
    txt = 'YWB1(J) = -YWD(J)'
    if LB > 0:
        txt += ' + YWB1(J)*(NE(J+1) - NE(J))/HRO'
    if LC > 0:
        txt += ' + YWC( J)*(TI(J+1) - TI(J))/HRO'
    teti += txt + '\n'
    if 'DVE' in var_defined:
        teti += 'YWB1(J) = YWB1(J) + DVE(J)*0.5*(NE(j) + NE(j+1))*log(TE(j)/TE(j+1))/HRO\n'
    teti += 'enddo\n'  # model1.f90, line 2292

    if 'TE' not in var_defined:
        logger.warning('Initial condition for TE is not defined\nTE=TEX(TSTART) will be used')

    teti += bnd_init('TE', var_defined, parse).replace('YWC', 'YWC1').replace('bctype', 'bc_type_imp(1)')

    if 'DSE' in var_defined:
        teti += 'DSE(ND1) = 1.\n'
    else:
        teti += 'DSE(ND1) = 0.\n'

    teti += const_text.TETIEQN.eqn

#---------------------------
# Ti part in implicit scheme
#---------------------------

    assign_type = parse.assign_d['TI']
    impl, asstyp = assign_type.split('_', 1)

    if none_in(['DI', 'HI', 'XI', 'CI', 'PI', 'PIT'], var_defined):
        for lbl in ('TIB', 'QIB', 'QITB'):
            if lbl in var_defined:
                logger.warning('Equation for "TI" requested but not defined')
                logger.warning('Boundary condition ignored')
                print(pa.apptmp(lbl, parse))
                var_defined.remove(lbl)

    if asstyp[:2] == 'EQ':
        teti += '! **** Ion temperature equation\n'
        teti += 'call markloc("TI equation")\n'
    else:
        logger.error('ERROR tetieqn: TI cannot be ASSIGNED with implicit scheme')

    teti += 'do j=1, NA1\n'
    for lbl in ['DI', 'HI', 'XI', 'CI', 'PI', 'PIT', 'DVI', 'DSI']:
        if lbl in var_defined:
            teti += pa.apptmp(lbl, parse)
    teti += 'if (j > NA) CYCLE\n'
    txt = 'YWA2(J) = 0.'
    if 'XI' in var_defined:
        txt += ' + XI(J)'
    if 'XN' in var_defined:
        txt += ' + GN2I*XN(J)'
    if 'DVI' in var_defined:
        txt += ' + DVI(J)'
    teti += txt + '\n'
    teti += 'YWA2(J) = YWA2(J)*(NI(J+1) + NI(J))*0.5\n'

    LB = 1
    if 'DN' in var_defined:
        if 'DI' in var_defined:
            teti += 'YWB2(J) = (DI(J) + GN2I*DN(J))\n'
        else:
            teti += 'YWB2(J) = GN2I*DN(J)\n'
    else:
        if 'DI' in var_defined:
            teti += 'YWB2(J) = DI(J)\n'
        else:
            LB = 0

    LC = 1
    if 'HN' in var_defined:
        if 'HI' in var_defined:
            teti += 'YWC(J) = (HI(J) + GN2I*HN(J))*(NI(J+1) + NI(J))/(TE(J+1) + TE(J))\n'
        else:
            teti += 'YWC(J) = GN2I*HN(J)*(NI(J+1) + NI(J))/(TE(J+1) + TE(J))\n'
    else:
        if 'HI' in var_defined:
            teti += 'YWC(J) = HI(J)*(NI(J+1) + NI(J))/(TE(J+1) + TE(J))\n'
        else:
            LC = 0

    txt = 'YWD(J) = 2.*GN2I*GNX(J)*SLAT(J)/G11(J)/(NE(J+1) + NE(J))'
    if 'CI' in var_defined:
        txt += ' + CI(J)'
    if 'CN' in var_defined:
        txt += '+ GN2I*CN(J)'
    teti += txt + '\n'
    teti += 'YWD(J) = YWD(J)*(NI(J+1) + NI(J))*0.5\n'

    teti += 'enddo\n'

    teti += 'do J=1, NA1\n'
    if 'PIT' not in var_defined:
        teti += 'PIT(J) = 0.\n'
    if 'PI' not in var_defined:
        teti += 'PI(J) = 0.\n'
    teti += 'PITOT(J) = PI(J)\n'
    teti += 'if (j > NA) CYCLE\n'
    txt = 'YWB2(J) = -YWD(J)'
    if LB > 0:
        txt += ' + YWB2(J)*(NE(J+1) - NE(J))/HRO'
    if LC > 0:
        txt += ' + YWC( J)*(TE(J+1) - TE(J))/HRO'
    teti += txt + '\n'
    if 'DVI' in var_defined:
        teti += 'YWB2(J) = YWB2(J) + DVI(J)*0.5*(NI(j) + NI(j+1))*log(TI(j)/TI(j+1))/HRO\n'
    teti += 'enddo\n'  # model1.f90, line 2547

    if 'TI' not in var_defined:
        logger.warning('Initial condition for TI is not defined')
        logger.warning('Using TI=TIX(TSTART)')

    teti += bnd_init('TI', var_defined, parse).replace('YWC', 'YWC2').replace('bctype', 'bc_type_imp(2)')

    if 'DSI' in var_defined:
        teti += 'DSI(ND1) = 1.\n'
    else:
        teti += 'DSI(ND1) = 0.\n'

    teti += const_text.TETIEQN.runeq

    return teti
