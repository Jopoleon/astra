#!/usr/bin/env python

import os, sys, logging, argparse
import equ_parser, code_gen
from parse_as import write_fortran
import config

fmt = logging.Formatter('%(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('as_parse')
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)


def astra_parser(f_equ):

    parse = equ_parser.EQU_PARSER(f_equ)
    return code_gen.CODE_GEN(parse)


def write_tmp(txt, dir_out=None, fortran='f90'):

    if dir_out is None:
        dir_out = '.'
    f90_l  = [ 'astra_out', 'setvar', 'inivar', 'detvar', 'detvar_init', 'ininam', 'subproc', 'postep', 'init_converge_step', 'eqns_inc']
    for lbl in f90_l:
        f_f90 = '%s/%s.f90' %(dir_out, lbl)
        write_fortran(f_f90, txt.__dict__[lbl], fortran=fortran)

    write_fortran('%s/model.txt'  %dir_out, txt.mtxt, fortran=fortran)
    write_fortran('%s/declar.fml' %dir_out, txt.fml , fortran=fortran)
    write_fortran('%s/declar.fnc' %dir_out, txt.fnc , fortran=fortran)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra parser')
    parser.add_argument('-awd', help='ASTRA home dirpath', required=False, default=os.getenv('AWD'))
    parser.add_argument('-equ', help='ASTRA equ filepath', required=False)
    args = parser.parse_args()

    f_equ = args.equ
    if f_equ is None:
        if args.awd is None:
            logger.error('Need either argument -awd or -equ')
            sys.exit()
        f_equ = '%s/tmp/model.tmp' %args.awd
    config.awd = args.awd
    config.fml_dir = '%s/fml'  %args.awd
    config.fnc_dir = '%s/fnc'  %args.awd
    txt = astra_parser(f_equ)
    write_tmp(txt, dir_out='./tmp')
