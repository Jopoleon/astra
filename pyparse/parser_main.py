#!/usr/bin/env python3

import os, sys, logging, argparse
import equ_parser, code_gen, exp_parser
from parse_as import write_fortran
import config

fmt = logging.Formatter('%(name)s | %(levelname)s: %(message)s', '%H:%M:%S')
hnd = logging.StreamHandler()
hnd.setFormatter(fmt)
logger = logging.getLogger('as_parse')
logger.addHandler(hnd)
logger.setLevel(logging.DEBUG)


def astra_parser(f_equ, f_exp):

    equ = equ_parser.EQU_PARSER(f_equ)
    exp = exp_parser.EXP_PARSER(f_exp)
    for var in equ.arname:
        if var not in exp.profiles['label']:
            logger.warning('X array %s used in equ, but missing in exp\n', var)
# Sanity check for x-input array size
    for jvar, var in enumerate(exp.profiles['label']):
        if len(exp.profiles['rho']) > 500:
            logger.warning('X array %s has nrho > 500', var)
    nr_x_max = max([len(x) for x in exp.profiles['rho']])
    print("NRX_MAX", nr_x_max)
    return code_gen.CODE_GEN(equ, nr_x_max)


def write_tmp(txt, dir_out=None):

    if dir_out is None:
        dir_out = '.'
    f90_l  = [ 'associate_pointers', 'astra_out', 'inivar', 'detvar', 'set_graph_names', 'ininam', 'postep', 'init_converge_step', 'eqns_inc']
    for lbl in f90_l:
        f_f90 = '%s/%s.f90' %(dir_out, lbl)
        write_fortran(f_f90, txt.__dict__[lbl])

    write_fortran('%s/model.txt'  %dir_out, txt.mtxt)
    write_fortran('%s/declar.fml' %dir_out, txt.fml )
    write_fortran('%s/declar.fnc' %dir_out, txt.fnc )


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra parser')
    parser.add_argument('-equ', help='ASTRA equ filepath', required=True)
    parser.add_argument('-exp', help='ASTRA exp filepath', required=True)

    args = parser.parse_args()

    txt = astra_parser(args.equ, args.exp)
    write_tmp(txt, dir_out='./src/tmp')
