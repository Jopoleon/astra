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
        init_txt += 'call initipc(NB1)\n'
        eqns_txt += const_text.SUBPROC.sbp_init
        init_txt += const_text.SUBPROC.sbp_init
        init_txt += 'call markloc("inikids")\n'
        init_txt += 'call inikids(NSBP, 64, LISTSB)\n'
    eqns_txt += 'call markloc("eqns.inc")\n'
    init_txt += 'call markloc("init.inc")\n'
    str_out = 'NITOT = NITOT + 1\n'
    eqns_txt += str_out
    init_txt += str_out

# If Fj, UPAR missing, do not fall back to any default, just skip

    ne_as = eqns.neeqn(parse, assign_type='Missing')
    init_txt += ne_as
    if parse.assign_d['NE'] == 'Missing':
        eqns_txt += ne_as

    te_as = eqns.teeqn(parse, assign_type='Missing')
    init_txt += te_as
    if parse.assign_d['TE'] == 'Missing':
        eqns_txt += te_as

    ti_as = eqns.tieqn(parse, assign_type='Missing')
    init_txt += ti_as
    if parse.assign_d['TI'] == 'Missing':
        eqns_txt += ti_as

    if 'CU' in parse.var_defined:
        cu_as = eqns.cuasn(parse, itype=-1)
    elif 'MU' in parse.var_defined:
        cu_as = eqns.cuasn(parse, itype=0)
    else:
        cu_as = eqns.cuas_uloop(parse)
    init_txt += cu_as
    if parse.assign_d['CU'] == 'Missing':
        eqns_txt += cu_as

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
            eqns_txt += eqns.fjeqn(parse, jf, assign_type=parse.assign_d[fj])

    if parse.assign_d['NE'] != 'Missing':
        eqns_txt += eqns.neeqn(parse, assign_type=parse.assign_d['NE'])
        if config.checkeqn:
            import iondens_ass
            eqns_txt += iondens_ass.NIAS.iondensassign

    if 'implicit' in parse.assign_d['TE'] or 'implicit' in parse.assign_d['TI']:
        eqns_txt += eqns.tetieqn(parse)
    else:
        if parse.assign_d['TE'] != 'Missing':
            eqns_txt += eqns.teeqn(parse, assign_type=parse.assign_d['TE'])
        if parse.assign_d['TI'] != 'Missing':
            eqns_txt += eqns.tieqn(parse, assign_type=parse.assign_d['TI'])

    if parse.assign_d['UPAR'] != 'Missing':
        eqns_txt += eqns.upeqn(parse, assign_type=parse.assign_d['UPAR'])

    if parse.assign_d['CU'][:2] == 'EQ':
        eqns_txt += eqns.cueqn(parse)
    else:
        if 'CU' in parse.var_defined:
            eqns_txt += eqns.cuasn(parse, -1)
        elif 'MU' in parse.var_defined:
            eqns_txt += eqns.cuasn(parse, 0)
        else:
             eqns_txt += eqns.cuas_uloop(parse)

# Closing statements

    init_txt += const_text.INIT.end
    eqns_txt += 'call markloc("eqns done")\n'

    return eqns_txt, init_txt
