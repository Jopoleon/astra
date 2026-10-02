#!/usr/bin/env python3
"""Convert param_files/sparams.dat written by ASTRA 8 (sbr/a2strahl.f90, newer STRAHL
format) into the format read by the STRAHL of the ASTRA 7 package. Idempotent.

Differences handled (all found by running the old reader on the A8 file):
- "cv finite vol=1, finite diff=0" + one integer  ->  "max_step max_diff" (1000 -1.,
  the value used by every example parameter file of the old package);
- integer columns that A8 writes as reals: source block column 3 (rate from file),
  divertor-puff block columns 1 and 4 (divertor puff, prompt redeposition).
"""
import re, sys

def to_int(tok):
    return str(int(round(float(tok.replace('D', 'E').replace('d', 'e')))))

def main(path):
    lines = open(path).read().split('\n')
    out, mode = [], None
    for ln in lines:
        s = ln.strip()
        if ln.startswith('cv'):
            mode = None
            if 'finite vol=1, finite diff=0' in ln:
                out.append('cv     max. iterations at fixed time   stop iteration if change below(%)')
                mode = 'finvol'
                continue
            if re.match(r'cv +r_source-r_lcfs', ln):
                mode = 'source'
            elif re.match(r'cv +divertor puff', ln):
                mode = 'puff'
            out.append(ln)
            continue
        if mode == 'finvol' and s:
            out.append('          1000   -1.')
            mode = None
            continue
        if mode in ('source', 'puff'):
            if not s:
                mode = None
                out.append(ln)
                continue
            tok = s.split()
            try:
                if mode == 'source' and len(tok) >= 3:
                    tok[2] = to_int(tok[2])
                elif mode == 'puff' and len(tok) >= 4:
                    tok[0] = to_int(tok[0]); tok[3] = to_int(tok[3])
            except ValueError:
                pass
            out.append('  ' + '   '.join(tok))
            continue
        out.append(ln)
    open(path, 'w').write('\n'.join(out))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'param_files/sparams.dat')
