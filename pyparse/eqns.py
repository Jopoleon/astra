import logging
import const_text, config
import parse_as as pa

logger = logging.getLogger('as_parse.eqns')
logger.setLevel(logging.INFO)

# Input:
#    var_defined: list of all variables assigned in the equ file
#             (left hand side of any '=' statement)
#             "if 'XN' in var_defined" corresponds to "IF(NBEG(LXN).GT.0)"
#
#    ass_type/ITYPE: None/-1, 'AS'/0, 'EQ'/1, 'EQ[j,rho]'/1, 'FU'/?


def none_in(list_in, var_list):

    flag = True
    for lbl in list_in:
        if lbl in var_list:
            flag = False
            break

    return flag


def write_equ_bnd(prof, rhob):

    lin_int  = 'if (ND1 < NA) then\n'
    lin_int += 'YA = (%s(NA1) - %s(ND1))/(ROC - %s)\n' %(prof, prof, rhob)
    lin_int += 'YB = (%s(ND1)*ROC - %s(NA1)*%s)/(ROC - %s)\n' %(prof, prof, rhob, rhob)
    lin_int += \
'''do j=ND1+1, NA'
lin_int += '%s(j) = YB + RHO(j)*YA
enddo
endif
''' %prof


def cuasn(parse, itype, neq=1):
    '''
itype:
   -1: CU profile prescribed
    0: MU profile prescribed
    '''

    cuas_txt = const_text.CUAS.header
 
    if 'MV' in parse.var_defined:
        cuas_txt += const_text.CUAS.mv1
        cuas_txt += pa.apptmp('MV', parse)
        cuas_txt += const_text.CUAS.mv2
    else:
        cuas_txt += 'FV(1: NA1) = 0.\n'

    cuas_txt += 'do j=1, NA1\n'
    for var in ('DC', 'HC', 'XC', 'CD', 'CC'):
        cuas_txt += pa.apptmp(var, parse)
    if none_in(['HC', 'DC', 'XC'], parse.var_defined):
        if 'CUBS' in parse.var_defined:
            cuas_txt += pa.apptmp('CUBS', parse)
        else:
            cuas_txt += 'CUBS(J) = 0.\n'
    else:
        cuas_txt += const_text.CUAS.cubs
        cuas_txt += cubs(parse.var_defined)
        cuas_txt += 'endif\n'

    if 'CD' in parse.var_defined:
        cuas_txt += 'YWA(J) = CUBS(J) + CD(J)\n'
    else:
        cuas_txt += 'YWA(J) = CUBS(J)\n'

    cuas_txt += 'enddo ! j (radial loop)\n'
    if neq == 0:
        cuas_txt += 'FP(1) = FV(1)\n'
    else:
        cuas_txt += 'FP(1) = FPO(1) + UPL(1)*TAU\n'

    if itype == 0: # MU
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


def cuas_uloop(parse, neq=1):
    '''current adjusted to match prescribed Uloop'''

    cuasu_txt = const_text.CUAS.header
 
    if 'MV' in parse.var_defined:
        cuasu_txt += const_text.CUAS.mv1
        cuasu_txt += pa.apptmp('MV', parse)
        cuasu_txt += const_text.CUAS.mv2
    else:
        cuasu_txt += 'FV(1: NA1) = 0.\n'

    cuasu_txt += const_text.CUAS.uloop_1

    cuasu_txt += 'do j=1, NA1\n'
    for var in ('DC', 'HC', 'XC', 'CD', 'CC'):
        cuasu_txt += pa.apptmp(var, parse)
    if none_in(['HC', 'DC', 'XC'], parse.var_defined):
        if 'CUBS' in parse.var_defined:
            cuasu_txt += pa.apptmp('CUBS', parse)
        else:
            cuasu_txt += 'CUBS(J) = 0.\n'
    else:
        cuasu_txt += const_text.CUAS.cubs
        cuasu_txt += cubs(parse.var_defined)
        cuasu_txt += 'endif\n'

    if 'CD' in parse.var_defined:
        cuasu_txt += 'YWA(J) = CUBS(J) + CD(J)\n'
    else:
        cuasu_txt += 'YWA(J) = CUBS(J)\n'

    cuasu_txt += \
'''YWB(J) = YF*CC(J)*YB/IPOL(J)**2
YWD(J) = YWA(J)*YD/(IPOL(J)**3 * G33(J))
'''
    cuasu_txt += 'enddo ! j (radial loop)\n'

    cuasu_txt += const_text.CUAS.uloop_2
    if 'MV' in parse.var_defined:
        cuasu_txt += \
'''FP(NA1) = FP(NA)*(IPL + FPO(NA1)*FP(NA1)) + YDF
CV(NA1) = CV(NA) + YWD(NA1) + ((CV(NA) + YWD(NA1)) - (CV(NA-1) + YWD(NA)))
'''
    else:
        cuasu_txt += 'FP(NA1) = FP(NA)*(IPL + FPO(NA1)*FP(NA1))\n'
    cuasu_txt += const_text.CUAS.uloop_3

    if 'MV' in parse.var_defined:
        cuasu_txt += 'CU(1: NA1) = CU(1: NA1) - CV(1: NA1)\n'
    if neq == 0:
        cuasu_txt += const_text.CUAS.beta
    cuasu_txt += 'enddo ! JCALL\n'

    return cuasu_txt


def cubs(var_defined):

    cubs_line = 'CUBS(J) = YA*(FP(J+1) - FP(J))*( 0.'
    if 'HC' in var_defined:
        cubs_line += ' + HC(J)*(TE(J+1) - TE(J))/(TE(J+1) + TE(J))'
    if 'XC' in var_defined:
        cubs_line += ' + XC(J)*(TI(J+1) - TI(J))/(TI(J+1) + TI(J))'
    if 'DC' in var_defined:
        cubs_line += ' + DC(J)*(NE(J+1) - NE(J))/(NE(J+1) + NE(J))'
    cubs_line += ' )\n'

    return cubs_line


def pre_eqn(parse, key, assign_type=None):

    right_hand = parse.right_hand_d
    if assign_type is None:
        assign_type = parse.assign_d[key]

    var_defined = parse.var_defined

    if none_in(config.coeff_d[key] + config.flux_d[key], var_defined):
        for var in config.bnd_d[key]:
            if var in var_defined:
                logger.warning('Equation for "%s" requested, but not defined', key)
                logger.warning('Boundary condition ignored:')
                logger.warning('%s = %s', var, right_hand[var])
                var_defined.remove(var)
        if assign_type == 'Missing':
            return None, var_defined

    pre_txt = ''
    if assign_type in ('AS', 'Missing'):
        pre_txt += '! **** %s assignment\n' %config.labels_d[key]
        pre_txt += 'call markloc("%s assignment")\n' %key
    if assign_type[:2] == 'EQ':
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
        pre_txt += '%sTOT(J) = %s(J)\n' %(flux, flux) # misplaced?
        pre_txt += 'enddo\n'
        pre_txt += 'YWC(1) = HRO\n'

    return pre_txt, var_defined


def neeqn(parse, assign_type=None):

    key = 'NE'
    if assign_type is None:
        assign_type = parse.assign_d[key]

    short = config.short_d[key]

    ne_txt, var_defined = pre_eqn(parse, key, assign_type=assign_type)
    if ne_txt is None:
        return ''

# Previous neeqn

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

# End previous neeqn

    if assign_type != 'Missing':
        if key not in var_defined:
            logger.warning('Initial condition for %s is not defined', key)
            logger.warning('  Using %s=%sX(TSTART)', key, key)

    rho_bnd = pa.set_rho(assign_type)
    if rho_bnd is None:
        ne_txt += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        ne_txt += 'RO%s = %s' %(short, l2f)
        ne_txt += 'ND1 = NODE(RO%s)\n' %short

    ne_txt += 'NA1%s = ND1\n' %short
    ne_txt += 'ND = ND1 - 1\n'

    if assign_type == 'AS':
# code_gen.CODE_GEN.ne_as

        if key in var_defined:
            ne_txt += 'if (ND1 < NA1) then\n'
            if 'NEB' in var_defined:
                ne_txt += 'do j=ND1, NA\n'
            else:
                ne_txt += 'do j=ND1, NA1\n'
            ne_txt += pa.apptmp(key, parse)
            ne_txt += 'enddo\n'
            ne_txt += 'endif\n'
            ne_txt += 'do j=1, ND1\n'
            ne_txt += pa.apptmp(key, parse)
            ne_txt += 'enddo\n'

        ne_txt += const_text.NEEQN.assigned
        if 'SNN' in var_defined:
            ne_txt += 'SNTOT(NA) = SNTOT(NA) + SNN(NA)*NE(NA)\n'
        if ('QNB' in var_defined) and ('QNNB' not in var_defined) and \
           ('NEB' not in var_defined) and rho_bnd is None:
            ne_txt += 'QN(NA1) = QNB\n'
        else:
            ne_txt += 'QN(NA1) = QN(NA)\n'
        ne_txt += 'GN(NA1) = QN(NA1)/SLAT(NA1)\n'
        ne_txt += 'SNTOT(NA1) = SNTOT(NA)\n'

    elif assign_type[:2] == 'EQ':      
# code_gen.CODE_GEN.ne_eq
        if key in var_defined:
            ne_txt += 'if (ND1 < NA1) then\n'
            ne_txt += 'do j=ND1, NA1\n'
            ne_txt += pa.apptmp(key, parse)
            ne_txt += 'enddo\n'
            ne_txt += 'endif\n'

        bnd_count = 0
        for var in config.bnd_d[key]:
            if var in var_defined:
                bnd_count += 1
        if bnd_count == 0: # No boundary condition is set
            logger.warning('Boundary condition for %s is not given', key)
            if rho_bnd is None:
                logger.warning('  It is set to %sX(t0)', key)
            elif key not in var_defined:
                logger.warning('  Using %sX(t0) at the shifted boundary', key)
            else:
                logger.warning('  Using %sX(t) at the shifted boundary', key)
            ne_txt += 'NEO(ND1: NA1) = NE(ND1: NA1)\n'
            ne_txt += 'YWC(4) = 1.\n'
        else:
            if 'NEB' in var_defined:
                ne_txt += 'NE(ND1) = %s' %pa.LINE2FOR(parse.right_hand_d['NEB'], parse)
                ne_txt += 'NEO(ND1: NA1) = NE(ND1: NA1)\n'
                ne_txt += 'YWC(4) = 1.\n'
            else:
                j_var = 2
                for var in ('QNB', 'QNNB'):
                    if var in var_defined:
                        ne_txt += pa.apptmp(var, parse)
                        ne_txt += 'YWC(%d) = %s\n' %(j_var, var)
                    else:
                        ne_txt += 'YWC(%d) = 0.\n' %j_var
                    j_var += 1
                ne_txt += 'YWC(4) = -1.\n'
    
        ne_txt += const_text.NEEQN.eqn

    ne_txt += 'do J=1, NA\n'
    if assign_type[:2] != 'EQ':
        ne_txt += 'QN(J) = G11(J)*(-YWA(J)*(NE(J+1) - NE(J))/HRO - 0.5*YWB(J)*(NE(J+1) + NE(J)))\n'

    ne_txt += 'QN(J) = QN(J) + SLAT(J)*GNX(J)\n'
    ne_txt += 'GN(J) = QN(J)/SLAT(J)\n'
    if 'SNN' in var_defined:
        ne_txt += 'SNTOT(J) = SNTOT(J) + SNN(J)*NE(J)\n'
    ne_txt += 'enddo\n'

    if ('QNB' in var_defined) and ('QNNB' not in var_defined) and \
       ('NEB' not in var_defined) and rho_bnd is None:
        ne_txt += 'QN(NA1) = QNB\n'
    else:
        ne_txt += 'QN(NA1) = QN(NA)\n'
    ne_txt += 'GN(NA1) = QN(NA1)/SLAT(NA1)\n'
    ne_txt += 'SNTOT(NA1) = SNTOT(NA)\n'

    return ne_txt


def tieqn(parse, assign_type=None):

    key = 'TI'
    if assign_type is None:
        assign_type = parse.assign_d[key]

    short = config.short_d[key]

# Previous tieqn

    ti_txt, var_defined = pre_eqn(parse, key, assign_type=assign_type)
    if ti_txt is None:
        return ''

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
        ti_txt += 'YWD(J) = DVI(j)\n'
    ti_txt += 'YWC(J) = PI(j)\n'
    ti_txt += 'enddo\n'

# End previous tieqn

    if (assign_type != 'Missing') and (key not in var_defined):
        logger.warning('Initial condition for %s is not defined', key)
        logger.warning('  Using %s=%sX(TSTART)', key, key)

    rho_bnd = pa.set_rho(assign_type)
    if rho_bnd is None:
        ti_txt += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        ti_txt += 'RO%s = %s' %(short, l2f)
        ti_txt += 'ND1 = NODE(RO%s)\n' %short

    ti_txt += 'NA1%s = ND1\n' %short
    ti_txt += 'ND = ND1 - 1\n'

    if assign_type == 'AS':
# code_gen.CODE_GEN.t_as

        if key in var_defined:
            ti_txt += 'if (ND1 < NA1) then\n'
            if 'TIB' in var_defined:
                ti_txt += 'do j=ND1, NA\n'
            else:
                ti_txt += 'do j=ND1, NA1\n'
            ti_txt += pa.apptmp(key, parse)
            ti_txt += 'enddo\n'
            ti_txt += 'endif\n'
            ti_txt += 'do j=1, ND1\n'
            ti_txt += pa.apptmp(key, parse)
            ti_txt += 'enddo\n'
    
        ti_txt += const_text.TIEQN.assigned
        if ('QIB' in var_defined) and ('QITB' not in var_defined) and \
           ('TIB' not in var_defined) and rho_bnd is None:
            ti_txt += 'QIB\n'
        else:
            ti_txt += 'QI(NA)\n'

    elif assign_type[:2] == 'EQ':
# code_gen.CODE_GEN.t_eq

        ti_txt += 'QI(1)  = RHO(ND1) - RHO(ND)\n'
        ti_txt += 'YWC(1) = RHO(ND1) - RHO(ND)\n'
        if key in var_defined:
            ti_txt += 'if (ND1 < NA1) then\n'
            if 'TIB' in var_defined:
                ti_txt += 'do j=ND1, NA\n'
            else:
                ti_txt += 'do j=ND1, NA1\n' # NA1, not NA
            ti_txt += pa.apptmp(key, parse)
            ti_txt += 'enddo\n'
            ti_txt += 'endif\n'

        bnd_count = 0
        for var in config.bnd_d[key]:
            if var in var_defined:
                bnd_count += 1
        if bnd_count == 0:
            logger.warning( 'Boundary condition for TI is not given')
            if rho_bnd is None:
                logger.warning('  It is set to TIX(t0)')
            elif key not in var_defined:
                logger.warning('  Using TIX(t0) at the shifted boundary')
            else:
                logger.warning('  Using TIX(t) at the shifted boundary')
            ti_txt += 'TIO(ND1: NA1) = TI(ND1: NA1)\n'
            ti_txt += 'YWC(4) = 1.\n'
        else:
            if 'TIB' in var_defined:
                ti_txt += 'TI(ND1) = %s' %pa.LINE2FOR(parse.right_hand_d['TIB'], parse)
                ti_txt += 'TIO(ND1: NA1) = TI(ND1: NA1)\n'
                ti_txt += 'QI(4)  = 1.\n'
                ti_txt += 'YWC(4) = 1.\n'
            else:
                j_var = 2
                for var in ('QIB', 'QITB'):
                    if var in var_defined:
                        ti_txt += pa.apptmp(var, parse)
                        ti_txt += 'QI(%d)  = %s\n' %(j_var, var)
                        ti_txt += 'YWC(%d) = %s\n' %(j_var, var)
                    else:
                        ti_txt += 'QI(%d)  = 0.\n' %j_var
                        ti_txt += 'YWC(%d) = 0.\n' %j_var
                    j_var += 1
                ti_txt += 'QI(4)  = -1.\n'
                ti_txt += 'YWC(4) = -1.\n'
    
        if 'DVI' in var_defined:
            ti_txt += 'YWD(ND1) = 1.\n'
        else:
            ti_txt += 'YWD(ND1) = 0.\n'
        if 'DSI' in var_defined:
            ti_txt += 'DSI(ND1) = 1.\n'
        else:
            ti_txt += 'DSI(ND1) = 0.\n'

        ti_txt += const_text.TIEQN.eqn

    return ti_txt


def teeqn(parse, assign_type=None):

    key = 'TE'
    if assign_type is None:
        assign_type = parse.assign_d[key]

    short = config.short_d[key]

    te_txt, var_defined = pre_eqn(parse, key, assign_type=assign_type)
    if te_txt is None:
        return ''

# Previous teeqn

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
        te_txt +=  'YWD(J) = DVE(j)\n'
    te_txt += 'YWC(J) = PE(j)\n'
    te_txt += 'enddo\n'

# End previous teeqn

    if (assign_type != 'Missing') and (key not in var_defined):
        logger.warning('Initial condition for %s is not defined', key)
        logger.warning('  Using %s=%sX(TSTART)', key, key)

    rho_bnd = pa.set_rho(assign_type)
    if rho_bnd is None:
        te_txt += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        te_txt += 'RO%s = %s' %(short, l2f)
        te_txt += 'ND1 = NODE(RO%s)\n' %short

    te_txt += 'NA1%s = ND1\n' %short
    te_txt += 'ND = ND1-1\n'

    if assign_type == 'AS':
# code_gen.CODE_GEN.te_as

        if key in var_defined:
            te_txt += 'if (ND1 < NA1) then\n'
            if 'TEB' in var_defined:
                te_txt += 'do j=ND1, NA\n'
            else:
                te_txt += 'do j=ND1, NA1\n'
            te_txt += pa.apptmp(key, parse)
            te_txt += 'enddo\n'
            te_txt += 'endif\n'
            te_txt += 'do j=1, ND1\n'
            te_txt += pa.apptmp(key, parse)
            te_txt += 'enddo\n'

        te_txt += const_text.TEEQN.assigned
        if ('QEB' in var_defined) and ('QETB' not in var_defined) and \
           ('TEB' not in var_defined) and rho_bnd is None:
            te_txt += 'QEB\n'
        else:
            te_txt += 'QE(NA)\n'

    elif assign_type[:2] == 'EQ':
# code_gen.CODE_GEN.t_eq

        te_txt += 'QE(1)  = RHO(ND1) - RHO(ND)\n'
        te_txt += 'YWC(1) = RHO(ND1) - RHO(ND)\n'
        if key in var_defined:
            te_txt += 'if (ND1 < NA1) then\n'
            if 'TEB' in var_defined:
                te_txt += 'do j=ND1, NA\n'
            else:
                te_txt += 'do j=ND1, NA1\n' # NA1, not NA
            te_txt += pa.apptmp(key, parse)
            te_txt += 'enddo\n'
            te_txt += 'endif\n'
    
        bnd_count = 0
        for var in config.bnd_d[key]:
            if var in var_defined:
                bnd_count += 1
        if bnd_count == 0:
            logger.warning( 'Boundary condition for TE is not given')
            if rho_bnd is None:
                logger.warning('  It is set to TEX(t0)')
            elif key not in var_defined:
                logger.warning('  Using TEX(t0) at the shifted boundary')
            else:
                logger.warning('  Using TEX(t) at the shifted boundary')
                te_txt += 'TEO(ND1: NA1) = TE(ND1: NA1)\n'
                te_txt += 'YWC(4) = 1.\n'
        else:
            if 'TEB' in var_defined:
                te_txt += 'TE(ND1) = %s' %pa.LINE2FOR(parse.right_hand_d['TEB'], parse)
                te_txt += 'TEO(ND1: NA1) = TE(ND1: NA1)\n'
                te_txt += 'QE(4)  = 1.\n'
                te_txt += 'YWC(4) = 1.\n'
            else:
                j_var = 2
                for var in ('QEB', 'QETB'):
                    if var in var_defined:
                        te_txt += pa.apptmp(var, parse)
                        te_txt += 'QE(%d)  = %s\n' %(j_var, var)
                        te_txt += 'YWC(%d) = %s\n' %(j_var, var)
                    else:
                        te_txt += 'QE(%d)  = 0.\n' %j_var
                        te_txt += 'YWC(%d) = 0.\n' %j_var
                    j_var += 1
                te_txt += 'QE(4)  = -1.\n'
                te_txt += 'YWC(4) = -1.\n'
    
        if 'DVE' in var_defined:
            te_txt += 'YWD(ND1) = 1.\n'
        else:
            te_txt += 'YWD(ND1) = 0.\n'
        if 'DSE' in var_defined:
            te_txt += 'DSE(ND1) = 1.\n'
        else:
            te_txt += 'DSE(ND1) = 0.\n'

        te_txt += const_text.TEEQN.eqn

    return te_txt


def cueqn(parse):

    cueq_txt = const_text.CUEQN.header
    if 'MV' in parse.var_defined:
        cueq_txt += const_text.CUEQN.mv

    for coeff in ('DC', 'HC', 'XC', 'CD', 'CC'):
        cueq_txt += pa.apptmp(coeff, parse)

    if none_in(['DC', 'HC', 'XC'], parse.var_defined):
        if 'CUBS' in parse.var_defined:
            cueq_txt += pa.apptmp('CUBS', parse)
        else:
            cueq_txt += 'CUBS(J) = 0.\n'
    else:
        cueq_txt += 'if (j == NA1) then\n'
        cueq_txt += 'CUBS(J)= max(0., (2*CUBS(NA) - CUBS(NA-1)))\n'
        cueq_txt += 'else\n'
        cueq_txt += 'CUBS(J) = YA*(FP(J+1) - FP(J))*(0.'
        if 'HC' in parse.var_defined:
            cueq_txt += ' + HC(J)*(TE(J+1) - TE(J))/(TE(J+1) + TE(J))'
        if 'XC' in parse.var_defined:
            cueq_txt += ' + XC(J)*(TI(J+1) - TI(J))/(TI(J+1) + TI(J))'
        if 'DC' in parse.var_defined:
            cueq_txt += ' + DC(J)*(NE(J+1) - NE(J))/(NE(J+1) + NE(J))'
        cueq_txt += ' )*0.5*(IPOL(J) + IPOL(J+1))\n'
        cueq_txt += 'endif\n'

    if 'CD' in parse.var_defined:
        cueq_txt += pa.apptmp('CD', parse)
        cueq_txt += 'YWD(J) = CUBS(J) + CD(J)\n'
    else:
        cueq_txt += 'YWD(J) = CUBS(J)\n'

    cueq_txt += 'YWD(J) = YWD(J)*YD/(IPOL(J)**3*G33(J))\n'
    cueq_txt += 'YWB(J) = 0.4*GP*CC(J)*YC/IPOL(J)**2\n'
    cueq_txt += 'enddo\n'

    if 'UEXT' not in parse.var_defined and 'LEXT' not in parse.var_defined and 'IPL' in parse.var_defined:
        cueq_txt += const_text.CUEQN.prescribed_ipl
    if 'UEXT' in parse.var_defined and 'LEXT' not in parse.var_defined:
        cueq_txt += const_text.CUEQN.prescribed_uloop
    if 'LEXT' in parse.var_defined:
        cueq_txt += const_text.CUEQN.circuit_eqn

    cueq_txt += const_text.CUEQN.eqn

    if 'LEXT' in parse.var_defined:
        cueq_txt += 'PSIFB = PSIEXT - PSPLEX*dfpdrbm12\n'

    cueq_txt += 'do J=1, NA1\n'

    if 'MV' in parse.var_defined:
        cueq_txt += 'CU(J) = CU(J) - CV(j)\n'

    cueq_txt += \
'''UPL(J) = (FP(J) - FPO(J))/TAU - YQDCMF(J)
ULON(J) = IPOL(J)*G33(J)*UPL(J)
UPL(J) = UPL(J) + YWR(1)
enddo
call CUOFP
'''

    if 'UEXT' in parse.var_defined or 'LEXT' in parse.var_defined:
        cueq_txt += '!DFPDR = ((FP(NA1) - FP(NA) - (FV(NA1) - FV(NA)))/HRO)\n'
        cueq_txt += 'IPL = 5.*IPOL(NA1)*G22(NA)*((FP(NA1) - FP(NA) - (FV(NA1) - FV(NA)))/HRO)/GP2/RTOR\n'

    return cueq_txt


def fjeqn(parse, jeq, assign_type=None):

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

    if assign_type is None:
        assign_type = parse.assign_d[key]

    short = config.short_d[key]

    fj_txt, var_defined = pre_eqn(parse, key, assign_type=assign_type)
    if fj_txt is None:
        return ''

# Previous fjeqn

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

# end previous fjeqn

    if assign_type != 'Missing':

        if  key not in var_defined:
            logger.warning('Initial condition for %s is not defined', key)
            logger.warning('  Using %s=%sX(TSTART)', key, key)

        rho_bnd = pa.set_rho(assign_type)
        if rho_bnd is None:
            fj_txt += 'ND1 = NA1\n'
        else:
            l2f = pa.LINE2FOR(rho_bnd, parse)
            fj_txt += '%s = %s' %(ro, l2f)
            fj_txt += 'ND1 = NODE(%s)\n' %ro

        fj_txt += 'NA1%s = ND1\n' %short
        fj_txt += 'ND = ND1 - 1\n'

        if assign_type[:2] == 'EQ':
            fj_txt += 'YWC(1) = HRO\n'

        if key in var_defined:
            fj_txt += 'if (ND1 < NA1) then\n'
            if assign_type == 'AS' and varb in var_defined:
                fj_txt += 'do j=ND1, NA\n'
            else:
                fj_txt += 'do j=ND1, NA1\n'
            fj_txt += pa.apptmp(key, parse)
            fj_txt += 'enddo\n'
            fj_txt += 'endif\n'

        if assign_type == 'AS':
            fj_txt += 'do j=1, ND1\n'
            fj_txt += pa.apptmp(key, parse)
            fj_txt += 'enddo\n'
        else:
            if qffb not in var_defined and qfb not in var_defined:
                if varb not in var_defined:
                    logger.warning('Boundary condition for %s is not set', key)
                else:
                    fj_txt += pa.apptmp(varb, parse)
                fj_txt += '%sO(ND1: NA1) = %s(ND1: NA1)\n' %(key, key)
                fj_txt += 'YWC(4) = 1.\n'

            if qfb in var_defined:
                fj_txt += pa.apptmp(qfb, parse)
                fj_txt += 'YWC(2) = %s\n' %qfb
                fj_txt += 'YWC(4) = -1.\n'
            else:
                fj_txt += 'YWC(2) = 0.\n'

            if qffb in var_defined:
                fj_txt += pa.apptmp(qffb, parse)
                fj_txt += 'YWC(3) = %s\n' %qffb
                fj_txt += 'YWC(4) = -1.\n'
            else:
                fj_txt += 'YWC(3) = 0.\n'

            fj_txt += const_text.FJEQN.eqn
            fj_txt += \
'call RUNEQ_EF(YWGN(1: NA1), YWHN(1: NA1), YWGO(1: NA1), YWHO(1: NA1), F%dO(1: NA1), YWNB(1: NA1), YWWB(1: NA1), YVR(1: NA1), YWM(1: NA1), G11(1: NA1), YWA(1: NA1), YWB(1: NA1), YWR(1: NA1), SFF%d(1: NA1),SF%d(1: NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ROC, RHO(1: NA1), imethod, YWC(1:7), F%d(1: NA1), QF%d(1: NA1), YQDCM(1: NA1), ADCMPF, MPHIT(1: NA1))\n' %(jeq, jeq, jeq, jeq, jeq)


    fj_txt += 'do j=1, NA-1\n'

    if assign_type in ('Missing', 'AS'):
        fj_txt += '%s(J) = -G11(J)*(YWA(J)*(%s(J+1) - %s(J))/HRO + 0.5*YWB(J)*(%s(J+1) + %s(J)))\n'      %(qf, key, key, key, key)
        if gf in var_defined:
            fj_txt += '%s(J) = %s(J) + SLAT(J)*%s(J)\n' %(qf, qf, gf)

    fj_txt += 'enddo\n'

    if gf not in var_defined:
        fj_txt += '%s(NA1) = %s(NA1)/SLAT(NA1)\n' %(gf, qf)

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


def upeqn(parse, assign_type=None):

    key = 'UPAR'
    if assign_type is None:
        assign_type = parse.assign_d[key]

    up_txt, var_defined = pre_eqn(parse, key, assign_type=assign_type)
    if up_txt is None:
        return ''

# Previous UPEQN

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
    up_txt += 'endif\n'

    up_txt += 'UPS0(j) = MRHO(j)*G41(j)*RTOR/IPOL(j)\n'
# neoclassical corrections
    up_txt += 'UPS1(j) = -DLNEO(j)*G41(j)*RTOR/IPOL(j)\n'
    up_txt += 'UPS2(j) = SGNEO(J)*RTOR*BTOR*IPOL(j)*(1. - BDB02(j)*G41(j)/IPOL(j)**0.2)\n'
    up_txt += 'enddo\n'

    if assign_type == 'Missing':
        up_txt += const_text.UPEQN.assigned
#        if 'TTRQB' in var_defined and 'UPARB' not in var_defined and 'ROU' not in var_defined:
#            up_txt += 'TTRQB\n'
#        else:
#            up_txt += 'QU(NA)\n'
        return up_txt

# From here, assign_type is either AS or EQ, not Missing

    rho_bnd = pa.set_rho(assign_type)
    if rho_bnd is None:
        up_txt += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        up_txt += 'ROU = %s' %l2f
        up_txt += 'ND1 = NODE(ROU)\n'

    if 'UPAR' not in var_defined:
        logger.warning('   Initial condition for UPAR is not defined\n' + \
                       '      UPAR=VTORX(TSTART) will be used')

    if assign_type == 'AS' and 'UPARB' in var_defined:
        up_txt += 'ND1 = NA1\n'
        up_txt += pa.apptmp('UPARB', parse)

#    if 'ROU' not in var_defined:
#        up_txt += 'ND1 = NA1\n'
#    else:
#        up_txt += pa.apptmp('ROU', parse)
#        up_txt += 'ND1 = NODE(ROU)\n'

    up_txt += 'NA1U = ND1\n'
    up_txt += 'ND = ND1 - 1\n'

    if assign_type[:2] == 'EQ':
        up_txt += 'QU(1)  = RHO(ND1) - RHO(ND)\n'
        up_txt += 'YWC(1) = RHO(ND1) - RHO(ND)\n'

    if 'UPAR' in var_defined:
        up_txt += 'if (ND1 < NA1) then\n'
        if assign_type == 'AS' and 'UPARB' in var_defined:
            up_txt += 'do j=ND1, NA\n'
        else:
            up_txt += 'do j=ND1, NA1\n'
        up_txt += pa.apptmp('UPAR', parse)
        up_txt += 'enddo\n'
        up_txt += 'endif\n'

    if assign_type == 'AS':
        up_txt += 'do j=1, ND1\n'
        up_txt += pa.apptmp('UPAR', parse)
        up_txt += 'enddo\n'
        up_txt += const_text.UPEQN.assigned
#        if 'TTRQB' in var_defined and 'UPARB' not in var_defined and 'ROU' not in var_defined:
#            up_txt += 'TTRQB\n'
#        else:
#            up_txt += 'QU(NA)\n'
        return up_txt

# From here, assign_type = EQ
# Line 1835 in model1.f90
    if 'TTRQB' not in var_defined and not 'UPARB' not in var_defined:
        logger.warning('Boundary condition for UPAR is not set')
        if 'ROU' not in var_defined:
            logger.warning('It is set to VTORX(t0)')
        else:
            if 'UPAR' in var_defined:
                logger.warning('VTORX(t) will be used at the shifted boundary')
            else:
                logger.warning('VTORX(t0) will be used at the shifted boundary')

    if 'UPARB' in var_defined:
        up_txt += pa.apptmp('UPARB', parse)
        if assign_type[:2] == 'EQ':
            up_txt += const_text.UPEQN.uparo
    elif 'TTRQB' in var_defined:
        up_txt += \
'''QU(2)  = TTRQB
YWC(2) = TTRQB
QU(3)  = 0.
QU(4)  = -1.
YWC(3) = 0.
YWC(4) = -1.
'''
    else:
        if assign_type[:2] == 'EQ':
            up_txt += const_text.UPEQN.uparo

    if assign_type[:2] == 'EQ':
        up_txt += const_text.UPEQN.eqn
    else:
        up_txt += const_text.UPEQN.assigned
        if 'TTRQB' in var_defined and 'UPARB' not in var_defined and 'ROU' not in var_defined:
            up_txt += 'QU(NA1) = TTRQB\n'
        else:
            up_txt += 'QU(NA1) = QU(NA)\n'

    return up_txt


def tetieqn(parse, itype=3):

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
    if 'HN' in var_defined and itype > 1:
        txt += ' + GN2E*HN(J)'
    if 'DVE' in var_defined:
        txt += ' + DVE(J)'
    teti += txt + '\n'
    teti += 'YWA1(J) = YWA1(J)*(NE(J+1) + NE(J))*0.5\n'

    LB = 1
    if itype > 1 and 'DN' in var_defined:
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
    if itype > 1 and 'XN' in var_defined:
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
    if 'CN' in var_defined and itype > 1:
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
        teti += 'YWD(J) = DVE(j)\n'
    teti += 'YWC(J) = PE(j)\n'
    teti += 'enddo\n'  # model1.f90, line 2292

    if 'TE' not in var_defined:
        logger.warning('Initial condition for TE is not defined\nTE=TEX(TSTART) will be used')

    rho_bnd = pa.set_rho(asstyp)
    if rho_bnd is None:
        teti += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        teti += 'ROE = %s' %(l2f)
        teti += 'ND1 = NODE(ROE)\n'
    teti += 'NA1E = ND1\n'
    teti += 'ND = ND1 - 1\n'

    teti += 'QE( 1) = RHO(ND1)-RHO(ND)\n'
    teti += 'YWC(1) = RHO(ND1)-RHO(ND)\n' # model1.f90, line 2315

    if 'TE' in var_defined:
        teti += 'if (ND1 < NA1) then\n'
        teti += 'do j=ND1, NA1\n'
        teti += pa.apptmp('TE', parse)
        teti += 'enddo\n'
        teti += 'endif\n'

# From here, itype is /= 0

    if none_in(['QETB', 'QEB', 'TEB'], var_defined):
        logger.warning('Boundary condition for TE is not set')
        if 'ROE' not in var_defined:
            logger.warning('It is set to TEX(t0)')
        else:
            if 'TE' in var_defined:
                logger.warning('Using TEX(t) at shifted boundary')
            else:
                logger.warning('Using TEX(t0) at shifted boundary')
        teti += const_text.TETIEQN.teold

    elif 'TEB' in var_defined:
        teti += 'TE(ND1) = %s' %pa.LINE2FOR(parse.right_hand_d['TEB'], parse)
        teti += const_text.TETIEQN.teold
    else:
        if 'QEB' in var_defined:
            teti += pa.apptmp('QEB', parse)
            teti += 'QE(2)   = QEB\n'
            teti += 'YWC1(2) = QEB\n'
        else:
            teti += 'QE(2)   = 0.\n'
            teti += 'YWC1(2) = 0.\n'
        if 'QETB' in var_defined:
            teti += pa.apptmp('QETB', parse)
            teti += 'QE(3)   = QETB\n'
            teti += 'YWC1(3) = QETB\n'
        else:
            teti += 'QE(3)   = 0.\n'
            teti += 'YWC1(3) = 0.\n'
        teti += 'QE(4)   = -1.\n'
        teti += 'YWC1(4) = -1.\n'

    if 'DVE' in var_defined:
        teti += 'YWD(ND1) = 1.\n'
    else:
        teti += 'YWD(ND1) = 0.\n'
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
    if 'XN' in var_defined and itype > 1:
        txt += ' + GN2I*XN(J)'
    if 'DVI' in var_defined:
        txt += ' + DVI(J)'
    teti += txt + '\n'
    teti += 'YWA2(J) = YWA2(J)*(NI(J+1) + NI(J))*0.5\n'

    LB = 1
    if itype > 1 and 'DN' in var_defined:
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
    if itype > 1 and 'HN' in var_defined:
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
    if 'CN' in var_defined and itype > 1:
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
        teti += 'YWD(J) = DVI(j)\n'
    teti += 'YWC(J) = PI(j)\n'
    teti += 'enddo\n'  # model1.f90, line 2547

    if 'TI' not in var_defined:
        logger.warning('Initial condition for TI is not defined')
        logger.warning('Using TI=TIX(TSTART)')

    rho_bnd = pa.set_rho(asstyp)
    if rho_bnd is None:
        teti += 'ND1 = NA1\n'
    else:
        l2f = pa.LINE2FOR(rho_bnd, parse)
        teti += 'ROI = %s' %l2f
        teti += 'ND1 = NODE(ROI)\n'
    teti += 'NA1I = ND1\n'
    teti += 'ND = ND1 - 1\n' # 2567
    teti += 'QI(1)  = RHO(ND1) - RHO(ND)\n'
    teti += 'YWC(1) = RHO(ND1) - RHO(ND)\n'
        
    if 'TI' in var_defined:
        teti += 'if (ND1 < NA1) then\n'
        teti += 'do j=ND1, NA1\n'
        teti += pa.apptmp('TI', parse)
        teti += 'enddo\n'
        teti += 'endif\n'

    if none_in(['QITB', 'QIB', 'TIB'], var_defined):
        logger.warning('Boundary condition for TI is not set')
        if 'ROI' in var_defined:
            if 'TI' in var_defined:
                logger.warning('Using TIX(t) at shifted boundary')
            else:
                logger.warning('Using TIX(t0) at shifted boundary')
        else:
            logger.warning('Using TIX(t0)')
        teti += const_text.TETIEQN.tiold
    elif 'TIB' in var_defined:
        teti += 'TI(ND1) = %s' %pa.LINE2FOR(parse.right_hand_d['TIB'], parse)
        teti += const_text.TETIEQN.tiold

    if 'QIB' not in var_defined and 'QITB' not in var_defined:
        teti += const_text.TETIEQN.tiold
    else:
        if 'QIB' in var_defined:
            teti += pa.apptmp('QIB', var_defined)
            teti += 'QI(2)   = QIB\n'
            teti += 'YWC2(2) = QIB\n'
        else:
            teti += 'QI(2)   = 0.\n'
            teti += 'YWC2(2) = 0.\n'
        if 'QITB' in var_defined:
            teti += pa.apptmp('QITB', var_defined)
            teti += 'QI(3)   = QITB\n'
            teti += 'YWC2(3) = QITB\n'
        else:
            teti += 'QI(3)   = 0.\n'
            teti += 'YWC2(3) = 0.\n'
        teti += 'QI(4)   = -1.\n'
        teti += 'YWC2(4) = -1.\n'

    if 'DVI' in var_defined:
        teti += 'YWD(ND1) = 1.\n'
    else:
        teti += 'YWD(ND1) = 0.\n'
    if 'DSI' in var_defined:
        teti += 'DSI(ND1) = 1.\n'
    else:
        teti += 'DSI(ND1) = 0.\n'

    teti += const_text.TETIEQN.runeq

    return teti
