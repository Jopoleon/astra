#!/usr/bin/env python3

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


def write_tmp(txt, dir_out=None):

    if dir_out is None:
        dir_out = '.'
    f90_l  = [ 'astra_out', 'setvar', 'inivar', 'detvar', 'detvar_init', 'ininam', 'init_sbp', 'shm2astra', 'postep', 'init_converge_step', 'eqns_inc']
    for lbl in f90_l:
        f_f90 = '%s/%s.f90' %(dir_out, lbl)
        write_fortran(f_f90, txt.__dict__[lbl])

    write_fortran('%s/model.txt'  %dir_out, txt.mtxt)
    write_fortran('%s/declar.fml' %dir_out, txt.fml )
    write_fortran('%s/declar.fnc' %dir_out, txt.fnc )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra parser')
    parser.add_argument('-equ', help='ASTRA equ filepath', required=True)

    args = parser.parse_args()

    txt = astra_parser(args.equ)
    write_tmp(txt, dir_out='./src/tmp')
