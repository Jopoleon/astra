#!/usr/bin/env python3

import os, sys, re, json, logging, argparse, time
import config, ufiles
import numpy as np

logger = logging.getLogger('exp_parser')
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)


def exp_split(exp):

    delimiters = ['POINTS', 'GRIDTYPE', 'FILTER', 'PROFIL', 'END']
    pattern = re.compile(r'\b(?:' + '|'.join(map(re.escape, delimiters)) + r')\b')

    match = pattern.search(exp)
    if match:
        exp1d = exp[:match.start()]
        exp2d = exp[match.start():]
    else:
        exp1d = exp
        exp2d = ''
    return exp1d, exp2d


def var_name(str_in):

    if not str_in or str_in[0] in (' ', '\t'): # If 1st character is tab or space, return blank
        return ''
    str_out = str_in.strip()
    if str_out[-1] == 'X':
        str_out = str_out[:-1]
    return str_out


def parse_u_line(str_in):

    var_name, str1 = str_in.split(' ', 1)
    words = str1.split(':')
    path = words[1].strip()
    if len(words) == 2:
        factor = 1.
    elif len(words) == 3:
        factor = float(words[2])
    return path, factor

    
class EXP_PARSER:


    def __init__(self, f_exp=None):

        tim = np.zeros(5)
        tim[0] = time.time()
        self.f_exp = f_exp

        f_json = '%s/astra_variables.json' %config.awd
        with open(f_json, 'r') as fjson:
            json_d = json.load(fjson)
        json_keys = {key: list(val.keys()) for key, val in json_d.items()}
        self.prof      = json_keys['profiles']
        self.profx     = json_keys['profiles_x']
        self.constants = json_keys['constants']
        self.intern1   = json_keys['internal']
        self.intern2   = json_keys['intern2']
        self.variables = json_keys['variables']
        self.varx      = json_keys['variables_x']

        self.scalars  = {}
        self.profiles = {}
        self.boundary = {}

        tim[1] = time.time()
        with open(f_exp, 'r') as f:
            exp = f.read()
        exp = re.sub(r'(?m)^\s*!.*\n?', '', exp) # Removes all lines starting with !
        exp = exp.upper()
        tim[2] = time.time()
        exp1d, exp2d = exp_split(exp)
        tim[3] = time.time()
        self.parse_exp1d(exp1d)
        tim[4] = time.time()
        print('Time analysis', np.diff(tim))


    def parse_exp1d(self, exp1d):

        lines = exp1d.splitlines()
        uf_keys = []

        for line in lines:
            varName = var_name(line[:6])
            tim   = line[8: 14]
            data  = line[16: 22]
            error = line[24: 30]
            if varName not in self.variables + ['NA1', 'TSTART', 'TEND']:
                continue
            if varName in ('NA1', 'TSTART', 'TEND', 'RTOR', 'AWALL', 'AB', 'ELONM', 'TRICH', 'AMJ', 'ZMJ'):
                if varName in self.scalars:
                    logger.error('>>> Error: faulty exp file %s: no time dependence allowed for variable %s', self.f_exp, varName)
                    sys.exit(1)
            if varName not in self.scalars:
                self.scalars[varName] = {'time': [], 'data': [], 'error': []}
            if varName in uf_keys: # double variable definition: 2x u-file, or u-file+exp_ascii
                logger.error('>>> Error: ambiguous definition of %s', varName)
                sys.exit(2)
            if 'U-FILE' in line: # u-file
                uf_path, factor = parse_u_line(line)
                f_in = '%s/udb/%s' %(config.awd, uf_path)
                print(f_in)
                uf = ufiles.UFILE(fin=f_in)
                self.scalars[varName]['time'] = uf.X['data']
                self.scalars[varName]['data'] = uf.f['data']
                uf_keys.append(varName)
            else: # ASCII block in exp file
                self.scalars[varName]['time'].append((float(tim) if tim.strip() else 0.))
                self.scalars[varName]['error'].append((float(error) if error.strip() else 0.))
                self.scalars[varName]['data'].append(float(data))
                

    def parse_exp2d(self, exp2d):

        lines = exp1d.splitlines()


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra parser')
    parser.add_argument('-exp', help='ASTRA exp filepath', required=False, default='%s/exp/aug34954_t' %config.awd)

    args = parser.parse_args()
    txt = EXP_PARSER(f_exp=args.exp)

    print(txt.scalars['NA1'])
    print(txt.scalars['AMJ'])
    print(txt.scalars['IPL'])
