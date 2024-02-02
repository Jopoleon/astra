import os, re, logging
import config

logger = logging.getLogger('as_parse.parse_as')
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)


def doublise(sarg):

    if '.' in sarg:
        if 'd' not in sarg and 'D' not in sarg and 'e' not in sarg and 'E' not in sarg:
            c = sarg.split('.')[0]
            sarg = sarg + 'd0'
        elif 'e' in sarg:
            sarg = sarg.replace('e', 'D')
        elif 'E' in sarg:
            sarg = sarg.replace('E', 'D')

    return sarg

#    try:
#        a = int(sarg)
#        return sarg
#    except:
#        try:
#            b = float(sarg)
#            c = sarg.split('.')[0]
#            return sarg + 'd0'
#        except:
#            return sarg


def apptmp(lbl, parse):

    pack = parse.fml_list, parse.fnc_list, parse.profiles, parse.arr_nam2
    txt = ''
    var = lbl.split('|', 1)[0]
    if lbl in parse.right_hand_d.keys():
        l2f = LINE2FOR(var, parse.right_hand_d[lbl], pack)
        txt += l2f.fcode

    return txt


def write_fortran(f_out, text, fortran='f77'):
    '''Writing FORTRAN tmp files'''

    indent_step = 4
    if fortran == 'f77':
        indent0 = 6
        line_break = '\n     &' + indent_step*' '
        llen = 68
    else:
        indent0 = 0
        line_break  = ' &\n' + indent_step*' '
        llen = 110 #90

    lines_in = text.split('\n')

# Add enddo-do in presence of "J+1" functions
    lines = []
    for line in lines_in:
        pieces = rec_split(line)
        if ('GRAD' in pieces) or ('GRADS' in pieces):
            lines.append('enddo')
            lines.append('do J=1, NA1')
        lines.append(line)

# Add proper indentation and line-breaking
    indent = indent0
    with open(f_out, 'w') as f:
        lin_old = ''
        for line in lines:
            lin_strip = line.lower().strip()
            lin_now = lin_strip.split('!')[0].strip()
            if len(lin_old.strip()) > 0:
                if lin_old[:2] == 'if' and lin_old[-4:] == 'then':
                    indent += indent_step
                if lin_old.split()[0].strip() == 'do':
                    indent += indent_step
                if lin_old == 'else':
                    indent += indent_step
            if lin_now in ('else', 'endif', 'enddo'):
                indent -= indent_step
            indent_str = indent*' '
            indented_line = indent_str + line
            if len(lin_strip) > 0:
                if lin_strip[0] == '!':
                   indented_line = line

            line_out = add_line_break(indented_line, llen=llen, line_break=line_break)
            lin_old = lin_now
            f.write(line_out)
            f.write('\n')
    logger.info('Written file %s' %f_out)


def add_line_break(line_in, llen=68, line_break='\n'):
    '''Splitting lines to be usable FORTRAN code'''

    if len(line_in) <= llen:
        line_out = line_in
    else:
        pieces = rec_split(line_in, syms='+ |- |* |/|(|)|,')
        line_out = ''
        leng_line = 0
        for piece in pieces:
            leng_line += len(piece)
            if leng_line <= llen:
                line_out += piece
            else:
                line_out += line_break + piece
                leng_line = len(line_break.split('\n')[1]) + len(piece)
    return line_out


def runeq(key, fl1, fl2, fl3):

    run_txt = """
call RUNEQ_EF(YWGN(1:NA1), YWHN(1:NA1), YWGO(1:NA1), YWHO(1:NA1), %sO(1:NA1), YWNB(1:NA1), YWWB(1:NA1), YVR(1:NA1), YWM(1:NA1), G11(1:NA1), YWA(1:NA1), """

    run_txt += 'YWB(1:NA1), YWR(1:NA1), %s(1:NA1), ' %fl1
    run_txt += '%s(1:NA1), RABDOT, BABDOT, ND1, NA1, HRO, TAU, ' %fl2
    run_txt += 'ROC, RHO(1:NA1), imethod, YWC(1:7), '
    run_txt += '%s(1:NA1), %s(1:NA1), YQDCM(1:NA1), ' %(key, fl3)
    run_txt += 'ADCMPF, MPHIT(1:NA1))\n'

    return run_txt


def set_rho(ass_type):

    rho_val = None
    if '[' in ass_type:
        tmp = ass_type.split('[')[1].split(']')[0]
        if ',' in tmp:
            tmp2 = tmp.split(',')
            jprof = int(tmp2[0])
            if jprof == 1:
                rho_val = 'RFAN(%s)' %tmp2[1]
            elif jprof == 2:
                rho_val = '%s*ROC' %tmp2[1]
            else:
                logger.error('First argument in ...:EQ[ , ] can be only either 1 or 2')
                logger.error('Please amend your equ file')
                sys.exit()
        else:
            rho_val = 'RFA(%s)' %tmp
    return rho_val


def insert_fml(fml):
    '''Replace "include fml/* with code inside'''

    txt = ''
    f_fml = '%s/%s' %(config.fml_dir, fml.lower())
    f = open(f_fml, 'r')
    for line in f.readlines():
        line = line.strip()
        if len(line) > 1:
            if line[0] != '!':
                txt += line.split('!')[0] + '\n'
    f.close()
    return txt.replace('INCLUDE', 'include')


def rec_split(line_in, syms='+|-|*|/|(|)|,'):
    '''Regex split with several delimiters'''

    syms2 = '('
    for sym in syms.split('|')[:-1]:
        syms2 += '\\' + sym + '|'
    syms2 += '\\' + syms.split('|')[-1] + ')'
    pieces = re.split(syms2, line_in)
    return [x for x in pieces if len(x.strip()) > 0]


def function_args(pieces):

    if pieces[0] != '(':
        return
    else:
        left_right_bra = -1
        for jpos, piece in enumerate(pieces[1:]):
            if piece == '(':
                left_right_bra -= 1
            if piece == ')':
                left_right_bra += 1
                if left_right_bra == 0:
                    break
        return jpos+2


def equ_prepare(f_equ):
    '''
    Skip blocks between "%" tags (comments)
    Skip lines beginning with "!" tag (comments)
    Trim lines
    Skip empty lines
    Put lines with ";" delimiter, such as "CAR1=CUECR; CAR2=PEECR;", into 2 lines
    '''

    try:
        fequ = open(f_equ, 'r', errors='replace')
    except:
        fequ = open(f_equ, 'r')

    equ_lines = []
    skip_block = False
    config.checkeqn = False

    for line in fequ.readlines():
        line = line.strip()
        if line == '':
            continue
        if line[0] == '%':
            skip_block = not skip_block
        elif not skip_block:
            if line[0] != '!': # Skip lines starting with '!'
                newline = line.split('!')[0] #Ignore text after '!' in a line
                for new_lin in newline.split(';'): # Split ';' into multiple lines
                    equ_lines.append(new_lin.strip())
            if line[:3] == '@!@':
                config.checkeqn = True

    fequ.close()
    return equ_lines


def indicise_lefteq(var, fnc_list, profiles, arr_nam2):
    '''Add proper FORTRAN index to ASTRA arrays, eqn left hand side'''

    var = var.strip()
    out = var
    tmp1 = var[:-1]
    if var in profiles + arr_nam2:
        out = '%s(J)' %var

    elif tmp1 in fnc_list + profiles:

        if var[-1] == 'B':
            if tmp1 in profiles:
                out = '%s(ND1)' %tmp1
        if var[-1] == 'C':
            out = '%s(1)' %tmp1

    return out


def indicise_righteq(line_in, fnc_list, profiles, arr_nam2):
    '''Add proper FORTRAN idnex to ASTRA arrays, eqn right hand side'''

    pieces = rec_split(line_in)

    line_out = ''
    for var in pieces:

        out = var      # including case var in ('+', '-', '*', '/', '(', ')', ',')
        if var in profiles + arr_nam2:
            out = '%s(J)' %var
        elif var in fnc_list:
            out = '%sR(RHO(J))' %var
        else:
            tmp1 = var[:-1]
            if var[-1] == 'B':
                if tmp1 in fnc_list + profiles:
                    if tmp1 in fnc_list:
                        out = '%sR(ROC)' %tmp1
                    else:
                        out = '%s(NA1)' %tmp1
            if var[-1] == 'C':
                if tmp1 in profiles:
                    out = 'RADIAL(%s, 0.d0)' %tmp1
                if tmp1 in fnc_list:
                    out = '%sR(0.d0)' %tmp1
        line_out += out
    return line_out


def format_number(str_num):
    '''
    Write numbers in double precision format
    If the input string is no number, it returns the input string unchanged
    '''

    try:
        num = float(str_num)
        str_out = '%e' %num
        while '0e' in str_out:
            str_out = str_out.replace('0e','e')
        str_out = str_out.replace('+0','')
        str_out = str_out.replace('-0','-')
        str_out = str_out.replace('e','d')
    except:
        str_out = str_num
    return str_out


def undef_inivar(var, short, defl='', varx2='', varx3=''):

    varx = var + 'X'
    if varx2 == '':
        varx2 = varx
    if varx3 == '':
        varx3 = varx + '(J)'
    inivar  = 'j1 = 0\n'
    inivar += 'do j=1, NARRX\n'
    inivar += 'if (EXARNM(j) == "%s" .and. IFDFAX(j) < 0) j1 = j\n' %varx2.ljust(6)
    inivar += 'enddo\n'
    inivar += 'if (NA1%s == NA1 .and. j1 /= 0 .and. ITREQ == 0) then\n' %short
    inivar += 'write(*, *) " >>> Warning: %s and %s are not defined"\n' %(var, varx)
    if defl != '':
        inivar += 'write(*, *) "         Using default value %s = %s"\n' %(var, defl)
    inivar += 'endif\n'
# Next lines: for all eqn_list variables
    inivar += 'if (IFDFAX(j1) > 0) then\n'
    inivar += 'do J=1, NA1\n'
    inivar += '%s(J) = %s\n' %(var, varx3)
    inivar += 'enddo\n'
    inivar += 'endif\n'

    return inivar


def write_declar(declar_files, ext='', ftype='double precision'):

    declar_txt = '%s :: ' %ftype

    for jvar, var in enumerate(declar_files):
        declar_txt += '%s%s' %(var, ext) # ext='R' for fnc, '' for fml
        if jvar < len(declar_files) - 1:
            declar_txt += ', '
    declar_txt += '\n'

    return declar_txt


def write_declar_fml(fml_files):

    dummy_flt = []
    dummy_int = []
    for fml in fml_files:
        ffml = '%s/%s' %(config.fml_dir, fml.lower())
        with open(ffml, 'r') as f:
            for lin in f.readlines():
                line = lin.strip().upper()
                trim = line.replace(' ', '')
                if trim == '':
                    continue
                if trim[0] == '!':
                    continue
 # some care required for lines with IF statements, that's why the split for ")"
                line1 = line
                if 'IF(' in trim:
                    line1 = line1.split(')')[1]
                if 'WHILE(' in trim:
                    line1 = line1.split(')')[1]
                if ('=' not in line1):
                    continue
                new_var = line1.split('=')[0].strip() # left-hand of a fortran statement
                if new_var.strip() == '':
                    continue
                if ' ' in new_var.strip():
                    new_var = new_var.split(' ')[-1].strip()
                if new_var.lower() in fml_files or new_var in fml_files:
                    continue
                if new_var in dummy_flt + dummy_int:
                    continue
                if new_var[0] in ('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z'):
                    dummy_flt.append(new_var)
                elif new_var[0] in ('I', 'J', 'K', 'L', 'M', 'N'):
                    dummy_int.append(new_var)
                else:
                    logger.warning('Variable name %s in file fml/%s is not allowed, skipped!' %(new_var, fml))

    logger.debug(dummy_flt)
    logger.debug(dummy_int)
    declar_txt  = write_declar(fml_files)
    if len(dummy_int) > 0:
        declar_txt += write_declar(dummy_int, ftype='integer')
    if len(dummy_flt) > 0:
        declar_txt += write_declar(dummy_flt)

    return declar_txt


def sbr_header(j_sbr, astr='', label='subroutine'):

    out_txt  = '! **** External %s\n' %label
    out_txt += 'IFSUB = 0\n'
    out_txt += 'if (KEY /= 0 .and. ABS(KEY - DTEQ(4,%d)) < 0.1) IFSUB = 1\n' %j_sbr
    out_txt += 'if (TIME >= DTEQ(2, %d) .and. TIME <= DTEQ(3, %d) .and. ' %(j_sbr, j_sbr)
    out_txt += 'TIME - TEQ(%d) + 1.E-7 - DTEQ(1, %d) > 0.0%s) IFSUB = 1\n' %(j_sbr, j_sbr, astr)

    return out_txt


def write_sbr(sbr_dic):

    out_txt = ''
    j_sbr  = sbr_dic['neq']
    sbrnam = sbr_dic['name']

    out_txt += 'if (IFSUB == 1) then\n'
    out_txt += 'TEQ(%d) = TIME\n' %j_sbr
    out_txt += 'call ADDTIME(cpuTime_tra)\n'
    out_txt += 'call markloc("subroutine %s")\n' %sbrnam
    if sbrnam == 'STRAHL':
        out_txt += 'call STRAHL(DTEQ(1, %d), %50s\n' %(j_sbr,'')
    else:
        out_txt += 'call %s(%s)\n' %(sbrnam, sbr_dic['args'])
    out_txt += 'call ADDTIME(cpuTime_sbr(%d))\n' %j_sbr
    out_txt += 'endif\n'

    return out_txt


def write_xpr(sbr_dic, j_ipc):

    out_txt = ''
    j_sbr  = sbr_dic['neq']
    sbrnam = sbr_dic['name']

    out_txt += 'if (IFSUB == 1) IFSUB = ifipc()\n'
    out_txt += 'if (IFSUB == 1) then\n'
    out_txt += 'TEQ(%d) = TIME\n' %j_sbr
    out_txt += 'IFSBP(%d) = %d\n' %(j_ipc, j_sbr)
    out_txt += 'call ADDTIME(cpuTime_tra)\n'
    sbr_nam = sbrnam.lower().replace('xpr/', '')
    out_txt += 'call markloc("subroutine %s")\n' %sbr_nam
    out_txt += 'call to_tra(%s, %d)\n' %(sbr_dic['args'], j_ipc)
    out_txt += 'call letsbp(%d)\n' %j_ipc
    out_txt += 'endif\n'

    return out_txt


def parse_inc(f_inc):

    profiles = []
    with open(f_inc, 'r') as f:
        for line in f.readlines()[2:]:
            arr = line.split()[0].strip()
            profiles.append(arr)

    return profiles


def parse_sbr(line):

# LOCSBR:
#   -4 Default value (the sbr is never used in the model)
#   -3 SB_P, tag "<" (detvar.tmp)
#   -2 SB_P, no tag  (init.inc & eqns.inc)
#   -1 SBR, tag "<"  (detvar.tmp)
#    0 SBR, no tag   (init.inc & eqns.inc)
#    1 SBR, tag ">"  (postep.inc)

    line = line.upper()
    tmp = line.replace(':=', '=').replace('=>', '=').strip()
    pieces = tmp.split(':')
    sbrnam = pieces[0].split('(')[0].strip('>').strip('<').strip()

    sbr_dic = {}
    args_str = ''
    if '(' in line:
        tmp1 = line.split('(', 1)[1]
        args_str = tmp1.rsplit(')', 1)[0]
    if line.count(':') == 4:
        dt, tmin, tmax, key = line.split(':')[1:]
    elif line.count(':') == 3:
        dt, tmin, tmax = line.split(':')[1:]
        key = ''
    elif line.count(':') == 2:
        dt, tmin = line.split(':')[1:]
        tmax = 1000.
        key = ''
    elif line.count(':') == 1:
        dt = line.split(':')[1]
        tmin = 0.
        tmax = 1000.
        key = ''
    if '&' in line: # SubProcess
        if '<' in line:
            locsbr = -3
        else:
            locsbr = -2
    else:           # subroutine
        if '<' in line:
            locsbr = -1  # detvar.f90; call in ASTRA_MAIN, STEPUP before equil
        elif '>' in line or sbrnam.upper() in ('MIXINT', 'MIXEXT', 'TSCTRL'):
            locsbr = 1   # postep.f90; call in STEPUP
        else:
            locsbr = 0   # init_converge_step.f90/eqns_inc.f90; call in ASTRA_MAIN, STEPUP

    args_str2 = ''
    if (args_str):
        args = args_str.split(',')
        for arg in args[:-1]:
            sarg = doublise(arg.strip())
            args_str2 += '%s, ' %sarg
        args_str2 += doublise(args[-1])

    sbr_dic['name'] = sbrnam
    sbr_dic['args'] = args_str2
    sbr_dic['tmin'] = tmin
    sbr_dic['tmax'] = tmax
    sbr_dic['dt']   = dt
    if key in ('', '>', '<'):
        sbr_dic['key']  = ''
    else:
        sbr_dic['key']  = ord(key.lower()) - ord("a") + 1
    sbr_dic['locsbr'] = locsbr

    return sbr_dic


class LINE2FOR:


    def __init__(self, left_hand, right_hand, pack):

        self.fml_list, self.fnc_list, self.profiles, self.arr_nam2 = pack
        self.flag_fml = {}
        for fml in self.fml_list:
            self.flag_fml[fml.lower()] = True
        if right_hand.strip() == '':
            self.fcode = '%s = 0.d0\n' %left_hand
        else:
            self.fcode = ''

            tmp = self.fmt_right_hand(right_hand)

            while 'replaced_fml' in tmp: # If there's more than one formula in one line
                self.fcode += tmp.replace('replaced_fml', '')
                tmp = self.fmt_right_hand(right_hand)
            right = tmp
            if left_hand.strip() == '':
                self.fcode = right
            else:
                left = indicise_lefteq(left_hand, self.fnc_list, self.profiles, self.arr_nam2)
                self.fcode += '%s = %s\n' %(left, right)


    def rec_fill_fml(self, txt_in):

        tmp1 = txt_in
        while('include' in tmp1):
            tout = ''
            lines = tmp1.split('\n')
            for line in lines:
                if 'include' in line:
                    fml = line.split('fml/')[1].split('\'')[0]

#                    if self.flag_fml[fml.lower()]:
#                        tout += insert_fml(fml)
#                        self.flag_fml[fml.lower()] = False
                    tout += insert_fml(fml)
                else:
                    if line.strip() != '':
                        tout += '%s\n' %line
            tmp1 = tout
        return tmp1.upper()


    def fmt_right_hand(self, line_in):
# Convert an "equ" statement in Fortran format

        tmp = line_in

        if line_in.strip() == '':
            line_out = ''
            return line_out

        logger.debug('IN:'+ line_in + '$')

# Check: exponential notation or real '-' between variables?
        for sym in ('E-', 'e-', 'D-', 'd-', 'E+', 'e+', 'D+', 'd+'):
            if sym[-1] == '-':
                repl = '$'
            else:
                repl = '#'
            if sym in tmp:
                piece0, piece1 = tmp.split(sym, 1)
                piece2 = rec_split(piece1)[0]
                piece3 = rec_split(piece0)[-1]
                logger.debug('piece3 %s', piece3)
                try:
                    flt = float(piece3) #Make sure there's a digit before "E-"
                    num = int(piece2) #Make sure there's an integer after "E-"
                    tmp = tmp.replace(sym, repl) # avoid splitting if exp notation
                except:
                    pass

        pieces = rec_split(tmp)

        line_out = ''
        n_pieces = len(pieces)

        jpos = 0
        while(jpos < n_pieces):

            var2 = pieces[jpos]
# avoid: integer array labels -> double precision
            var = format_number(var2).upper().strip()
            logger.debug('var2=%s, var=%s, jpos=%d', var2, var, jpos)
            if '.' not in var2:
                if pieces[jpos-1] == '(' and pieces[jpos+1] == ')' or \
                   pieces[jpos-1] == '(' and pieces[jpos+1] == ',' or \
                   pieces[jpos-1] == ',' and pieces[jpos+1] == ')':
                    var = var2

            if (jpos < n_pieces-2 and pieces[jpos + 2] == 'AFX'):
                count_close_bracket = 0
                block_left = ''
                jcomma = -1
                for jpiec, piec in enumerate(pieces[jpos+2:]):
                    if piec == '(':
                        count_close_bracket -= 1
                    elif piec == ')':
                        count_close_bracket += 1
                    elif piec == ',' and count_close_bracket == 0:
                        jcomma = jpiec
                    if count_close_bracket == 1:
                        break
                    if jcomma == -1:
                        block_left += format_number(piec)
                if block_left in self.fnc_list:
                    out = 'RADIAL(%sR, RFA(%s))' %(pieces[jpos], block_left)
                else:
                    out = 'RADIAL(%s, RFA(%s))' %(pieces[jpos], block_left)
                jpos += jpiec + 2

            elif var in ('VINT', 'IINT', 'LININT'):
                var3 = pieces[jpos+2]
                tmp3 = var3[:-1]
                if var3[-1] == 'B':
                    if tmp3 in self.fnc_list:
                        out = '%s(%sR, ROC)' %(var, tmp3)
                    if tmp3 in self.profiles:
                        out = '%s(%s, ROC)'  %(var, tmp3)
                    jpos += 3
                elif len(pieces[jpos+2:]) == 2:
                    logger.debug('var3 %s', var3)
                    if var3 in self.fnc_list:
                        out = '%s(%sR, j*HRO)' %(var, var3)
                    if var3 in self.profiles:
                        out = '%s(%s, j*HRO)'  %(var, var3)
                    jpos += 3
                else:
                    count_close_bracket = 0
                    block_left = ''
                    block_right = ''
                    jcomma = -1
                    for jpiec, piec in enumerate(pieces[jpos+2:]):
                        if piec == '(':
                            count_close_bracket -= 1
                        elif piec == ')':
                            count_close_bracket += 1
                        elif piec == ',' and count_close_bracket == 0:
                            jcomma = jpiec
                        if count_close_bracket == 1:
                            break
                        if jcomma == -1:
                            block_left += piec
                        elif jcomma != jpiec:
                            block_right += piec
                    if not block_right: # VINT(CAR11) * ...
                        block_right = 'j'
                    if block_left in self.fnc_list:
                        out = '%s(%sR, %s*ROC)' %(var, block_left, block_right)
                    else:
                        out = '%s(%s, %s*ROC)'  %(var, block_left, block_right)
                    jpos += jpiec + 2
                    logger.debug(block_left)
                    logger.debug(block_right)
                    logger.debug(out)

            elif var in ('RFVAL', 'RFVEX', 'RFVIN',
                         'AFVAL', 'AFVEX', 'AFVIN',
                         'ASTEP', 'RSTEP', 'XSTEP',
                         'ATX', 'ATR', 'V_95_POS', 'RFMAX', 'RFMIN', 'FRMAX', 'FRMIN'):
                jbra = function_args(pieces[jpos+1: n_pieces])
                out = ''
                for j in range(jpos, jpos+jbra):
                    if pieces[j] in self.profiles:
                        out += '%s(1:NA1)' %pieces[j]
                    else:
                        out += doublise(pieces[j])
                if var in ('RFVAL', 'RFVEX', 'RFVIN'):
                    out += '*ROC'
                elif var in ('AFVAL', 'AFVEX', 'AFVIN'):
                    out += '*ABC'
                elif var in ('ASTEP', 'RSTEP', 'XSTEP'):
                    out += ', J'
                out += ')'
                jpos += jbra

            elif var in ('GRAD', 'GRADS'):
                var3 = pieces[jpos+2]
                tmp3 = var3[:-1]
                if var3[-1] == 'B':
                    if tmp3 in self.fnc_list:
                        out = '%s(%sR, NA1)' %(var, tmp3)
                    if tmp3 in self.profiles:
                        out = '%s(%s, NA1)'  %(var, tmp3)
                elif var3[-1] == 'C':
                    if tmp3 in self.fnc_list:
                        out = '%s(%sR, 1)' %(var, tmp3)
                    if tmp3 in self.profiles:
                        out = '%s(%s, 1)'  %(var, tmp3)
                else:
                    if var3 in self.fnc_list:
                        out = '%s(%sR, J)' %(var, var3)
                    if var3 in self.profiles:
                        out = '%s(%s, J)'  %(var, var3)
                jpos += 3

# Formula
            elif var in self.fml_list:
                if self.flag_fml[var.lower()]:
                    tmp4 = insert_fml(var)
                    line_out = 'replaced_fml' + self.rec_fill_fml(tmp4)
                    self.flag_fml[var.lower()] = False
                    return line_out
                out = var
# Between brackets
            elif var in self.profiles + self.fnc_list:
                if (jpos < n_pieces - 3) and (pieces[jpos+1] == '(' and pieces[jpos+3] == ')'):
                    try: # double precision
                        arg = float(pieces[jpos+2])
                        out = 'RADIAL(%s, RFA(%s)' %(var, format_number(arg))
                        jpos += 2
                    except: # integer
                        out = var
                elif jpos < n_pieces - 5 and pieces[jpos+1] == '(' and pieces[jpos+3] != ')' and pieces[jpos+5] == ')':
                    try: # double precision
                        out = 'RADIAL(%s, RFA(%s%s%s)' %(var, format_number(pieces[jpos+2]), pieces[jpos+3], pieces[jpos+4])
                        jpos += 4
                    except: # integer
                        out = var
                else:
                    out = indicise_righteq(var, self.fnc_list, self.profiles, self.arr_nam2)
            else:
                out = indicise_righteq(var, self.fnc_list, self.profiles, self.arr_nam2)

            line_out += out
            jpos += 1

# Reinserting exponential notation, after parsing for operational '+', '-'
        line_out = line_out.replace('"', '')
        line_out = line_out.replace('$', 'd-')
        line_out = line_out.replace('#', 'd+')
        logger.debug('OUT: %s', line_out)
        logger.debug('')

        return line_out
