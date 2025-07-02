import const_text, eqns, config
from parse_as import *


def eqns_init(parse):
# Corresponds to model2:alldef

    eqns_lin = parse.eqns_lines

    eqns_txt = ''
    init_txt = ''

    sbp_lines = []
    sbr_txt = ''
    j_sbr = 1

    for line in parse.sbr_lines:

# Subroutines
        sbr_d = parse_sbr(line)
        sbr_d['neq'] = j_sbr
        if sbr_d['locsbr'] == 0: #sbr
            a_str = ''
#            if (sbr_d['name'] in ('MIXINT', 'MIXEXT', 'TSCTRL')):
#                a_str = '.and. JIT == JEX'
            sbr_txt += sbr_header(sbr_d['neq'], astr=a_str)
            sbr_txt += write_sbr(sbr_d)
# Subprocess
        if sbr_d['locsbr'] in (-2, -3):
            sbp_lines.append(line)
        j_sbr += 1

    NSBP = len(sbp_lines)
    NSBR = len(parse.sbr_lines)
    if NSBP == 0 and NSBR == 0:
        str_out = '! **** No external subroutines'
        eqns_txt += str_out
        init_txt += str_out
    if NSBP > 0:
        init_txt += 'call markloc("initipc")\n'
        init_txt += 'call initipc(NA1, n_ql)\n'
        eqns_txt += const_text.SUBPROC.sbp_init
        init_txt += const_text.SUBPROC.sbp_init
        init_txt += 'call markloc("inikids")\n'
        init_txt += 'call inikids(NSBP, 64, LISTSB)\n'
    eqns_txt += 'call markloc("eqns")\n'
    init_txt += 'call markloc("init")\n'
    str_out = 'NITOT = NITOT + 1\n'
    eqns_txt += str_out
    init_txt += str_out

# If Fj, UPAR missing, do not fall back to any default, just skip
    for lbl in ('NE', 'TE', 'TI'):
        lbl_init, _ = eqns.pre_eqn(parse, lbl, assign_type='Missing')
        init_txt += lbl_init

    if 'CU' in parse.var_defined:
        cu_as = eqns.cuasn(parse, bc='CU')
    elif 'MU' in parse.var_defined:
        cu_as = eqns.cuasn(parse, bc='MU')
    else:
        cu_as = eqns.cuas_uloop(parse)
    init_txt += cu_as

# Subroutines

    init_txt += sbr_txt
    eqns_txt += sbr_txt

# Subprocesses

    if NSBP > 0:
        statement = '\ncall SUBPROC\n\n'
        init_txt += statement
        eqns_txt += statement

# Equations

    for jf in range(10):
        fj = 'F%d' %jf
        if parse.assign_d[fj] != 'Missing':
            eqns_txt += eqns.fjeqn(parse, jf)

    eqns_txt += eqns.neeqn(parse)
    if config.checkeqn:
        import iondens_ass
        eqns_txt += iondens_ass.NIAS.iondensassign

    if 'implicit' in parse.assign_d['TE'] or 'implicit' in parse.assign_d['TI']:
        eqns_txt += eqns.tetieqn(parse)
    else:
        eqns_txt += eqns.teeqn(parse)
        eqns_txt += eqns.tieqn(parse)

    if parse.assign_d['UPAR'] != 'Missing':
        eqns_txt += eqns.upeqn(parse)

    if parse.assign_d['CU'][:2] == 'EQ':
        eqns_txt += eqns.cueqn(parse)
    else:
        eqns_txt += cu_as

    eqns_txt += const_text.pol_flux
    init_txt += const_text.pol_flux

# Closing statements

    init_txt += const_text.INIT.end
    eqns_txt += 'call markloc("eqns done")\n'

    return eqns_txt, init_txt
