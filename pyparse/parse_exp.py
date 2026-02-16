#!/usr/bin/env python3

import os, sys, re, json, logging, argparse, time
import config, ufiles
import numpy as np

logger = logging.getLogger('exp_parser')
#logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)


prof_attr_types = {
    'POINTS': int,
    'GRIDTYPE': int,
    'NTIMES': int,
    'FILTER': float,
    'FACTOR': float,
    'NAMEXP': str,
}
delimiters = list(prof_attr_types) + ['PROFILE', 'END']
grid_d = {'MINORRAD': 10, 'MAJORRAD': 19, 'RHO': 12, 'RHO_TOR': 12, 'POLOIDALFLUX': 13}


def exp_split(exp):
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
    return str_out.upper()


def append_x(str_in):
    if not str_in or str_in[0] in (' ', '\t'): # If 1st character is tab or space, return blank
        return ''
    else:
        return var_name(str_in) + 'X'


def parse_u_line(str_in):
    var_name, str1 = str_in.split(' ', 1)
    words = str1.split(':')
    path = words[1].strip()
    if len(words) == 2:
        factor = 1.
    elif len(words) == 3:
        factor = float(words[2])
    return var_name, path, factor


def parse_line2d(line):
    words = line.split()
    result = {}
    for i, word in enumerate(words[:-1]):   # avoid overflow
        if word in prof_attr_types:
            try:
                result[word] = prof_attr_types[word](words[i + 1])
            except ValueError:
                raise ValueError(f"Invalid value for {word}: {words[i+1]}")
    if 'NAMEXP' in result:
        result['NAMEXP'] = append_x(result['NAMEXP'])
    return result


def read_float_block(lines, start_line, end_line):
    block_lines = lines[start_line:end_line]   # end exclusive
    block_text = " ".join(block_lines)         # normalize spacing
    values = np.fromstring(block_text, sep=' ')
    return values


class EXP_PARSER:


    def __init__(self, f_exp=None):

        tim = np.zeros(7)
        tim[0] = time.time()
        self.f_exp = f_exp

        f_json = '%s/astra_variables.json' %config.awd
        with open(f_json, 'r') as fjson:
            json_d = json.load(fjson)
        json_keys = {key: list(val) for key, val in json_d.items()}
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
        self.coils    = {}

        tim[1] = time.time()
        with open(f_exp, 'r') as f:
            exp = f.read()
        exp = re.sub(r'(?m)^\s*!.*\n?', '', exp) # Removes all lines starting with !
        tim[2] = time.time()
        exp1d, exp2d = exp_split(exp)
        tim[3] = time.time()
        self.parse_exp1d(exp1d)
        tim[4] = time.time()
        self.parse_exp2d(exp2d)
        tim[5] = time.time()
        self.write_json()
        tim[6] = time.time()
        print('Time analysis', np.diff(tim))


    def parse_exp1d(self, exp1d):
# Input is already uppercase

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
            if 'U-FILE' in line.upper(): # u-file
                varName, uf_path, factor = parse_u_line(line)
                f_in = '%s/udb/%s' %(config.awd, uf_path)
                uf = ufiles.UFILE(fin=f_in)
                self.scalars[varName]['time'] = uf.X['data'].tolist()
                self.scalars[varName]['data'] = (factor*uf.f['data']).tolist()
                self.scalars[varName]['error'] = len(self.scalars[varName]['data'])*[0.]
                uf_keys.append(varName)
            else: # ASCII blocks in exp file
                self.scalars[varName]['time'].append((float(tim) if tim.strip() else 0.))
                self.scalars[varName]['error'].append((float(error) if error.strip() else 0.))
                self.scalars[varName]['data'].append(float(data))
                

    def parse_exp2d(self, exp2d):
# Input is already uppercase

        varInExp = [append_x(x) for x in re.findall(r'\bNAMEXP\s+(\S+)', exp2d)]
        for var in ('BNDUX', 'BNDX', 'CCOILX', 'VCOILX'):
            if varInExp.count(var) > 1:
                logger.error('Variable %s is defined more than once in the exp file', var)
                sys.exit(3)

        lines = exp2d.splitlines()
        n_lines = len(lines)
        keywords = delimiters + ['U-file']
        pattern = re.compile(r'\b(?:' + '|'.join(keywords) + r')\b')
        line_numbers = [i for i, line in enumerate(lines) if pattern.search(line)]

        uf_keys = []
        self.boundary['time'] = []
        for lbl in ('label', 'grid_type', 'filter', 'time', 'rho', 'data'): 
             self.profiles[lbl] = []

        alpha_glob = 0.001

        for j, jlin in enumerate(line_numbers): # jumping between "string" lines, either exp-ASCII blocks or u-file lines
            line = lines[jlin]
            line_strip = line.strip()
            if not line_strip:
                continue
            if line_strip.startswith('!'):
                continue
            if line_strip.startswith('END'):
                break

            alpha = alpha_glob
            print('LINE:', line)
            if j < len(line_numbers) - 1:
                jnext = line_numbers[j+1]
            else:
                jnext = n_lines

            if 'U-file' in line_strip: # u-file
                varName, uf_path, factor = parse_u_line(line)
                f_in = '%s/udb/%s' %(config.awd, uf_path)
                uf = ufiles.UFILE(fin=f_in)
                gridtype = grid_d[uf.Y['label'].strip().upper()]
                varName = append_x(varName)
                self.profiles['grid_type'].append(gridtype)
                self.profiles['label'].append(varName)
                self.profiles['time'].append(uf.X['data'].tolist())
                self.profiles['rho'].append(uf.Y['data'].tolist())
                self.profiles['data'].append((factor*uf.f['data']).ravel().tolist())
                self.profiles['filter'].append(alpha)
                uf_keys.append(varName)
                
            else: # A line of the kind "NAMEXP ...  GRIDTYPE ...", not u-file
                attr_d = parse_line2d(line)
                if 'FILTER' in attr_d:
                    alpha = attr_d['FILTER']
                    if line.startswith('FILTER'):
                        alpha_glob = alpha
                if 'NAMEXP' not in attr_d:
                    continue

                varName = attr_d['NAMEXP'] # already with trailing 'X'
                dataStream = read_float_block(lines, jlin+1, jnext)
                if varName in ('CCOILX', 'VCOILX'): # differs from (say) NEX because it has no radial grid
                    if 'NTIMES' in attr_d:
                        nt = attr_d['NTIMES']
                    else:
                        nt = 1
                    self.coils[varName.lower()] = {
                        'time': dataStream[:nt].tolist(),
                        'data': dataStream[nt:].tolist()}

                elif varName == 'BNDX':
                    if self.boundary['time']:
                        logger.error('Boundary must be defined in a single group')
                        sys.exit(3)
                    if 'POINTS' in attr_d:
                        nthe = attr_d['POINTS']
                    else:
                        logger.error('Number of boundary points must be defined')
                        sys.exit(4)
                    if 'NTIMES' in attr_d:
                        nt = attr_d['NTIMES']
                    else:
                        nt = 1
                    self.boundary['time'] = dataStream[:nt].tolist()
                    bnd_rz = dataStream[nt:]
                    assert bnd_rz.size == 2 * nt * nthe
                    tmp = bnd_rz.reshape(nthe, 2, nt)
                    self.boundary['R'] = tmp[:, 0, :].ravel().tolist()
                    self.boundary['Z'] = tmp[:, 1, :].ravel().tolist()
                elif varName == 'BNDUX':
                    uf_path = lines[jlin+1].strip()
                    fR_in = '%s/udb/%s_r' %(config.awd, uf_path)
                    fZ_in = '%s/udb/%s_z' %(config.awd, uf_path)
                    ufR = ufiles.UFILE(fin=fR_in)
                    ufZ = ufiles.UFILE(fin=fZ_in)
                    self.boundary['time'] = ufR.X['data'].tolist()
                    self.boundary['R'] = ufR.f['data'].ravel().tolist()
                    self.boundary['Z'] = ufZ.f['data'].ravel().tolist()

                elif varName in self.profx:
                    if 'NTIMES' in attr_d:
                        nt = attr_d['NTIMES']
                    else:
                        nt = 1
                    nrho = attr_d['POINTS']
                    self.profiles['label'].append(varName)
                    self.profiles['grid_type'].append([attr_d['GRIDTYPE']])
                    if 'FILTER' in attr_d:
                        self.profiles['filter'].append([attr_d['FILTER']])
                    else:
                        self.profiles['filter'].append(nt*[alpha_glob])
                    if attr_d['GRIDTYPE'] < 17:
                        assert dataStream.size == nt + nrho + nt*nrho
                        time, rho, data = np.split(dataStream, [nt, nt+nrho])
                    else:
                        assert dataStream.size == nt + 2*nrho + nt*nrho
                        time, rho, data = np.split(dataStream, [nt, nt+2*nrho])
                    self.profiles['time'].append(time.tolist())
                    self.profiles['rho'].append(rho.tolist())
                    self.profiles['data'].append(data.tolist())


    def write_json(self):

        data = {"coils": self.coils, "scalars": self.scalars, "boundary": self.boundary, "profiles": self.profiles}
        with open('data.json', 'w') as f:
            json.dump(data, f, indent=2)

        
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='astra parser')
    parser.add_argument('-exp', help='ASTRA exp filepath', required=False, default='%s/exp/aug34954_t' %config.awd)

    args = parser.parse_args()

#    args.exp='%s/exp/30000_3.4' %config.awd

    raw = EXP_PARSER(f_exp=args.exp)

    if 'ccoilx' in raw.coils:
        print('CCOILX')
        print(raw.coils['ccoilx']['time'])
        print(raw.coils['ccoilx']['data'])
    print(args.exp)
    print(raw.profiles['label'])
    n_len = len([x for xs in raw.profiles['data'] for x in xs])
    print('Size of profiles array', n_len)
