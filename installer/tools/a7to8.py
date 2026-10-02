#!/usr/bin/env python3
"""Bring an ASTRA 7 user area (models, formulas, data) into an ASTRA 8 installation.

    a7to8.py <A7 user dir> <A8 dir> MODEL [MODEL ...] [--exp FILE ...] [--dry-run]

For every MODEL (file name in <A7>/equ):
  * equ/<name>      model copy; quoted A7 Fortran "work(j,N)" written by NEOCL4 is replaced
                    by the ASTRA 8 nclass_mod arrays (A8 NEOCL4 no longer fills WORK);
  * fml/<f>         formulas used by the model (and their includes) that ASTRA 8 does not have,
                    converted from A7 fixed form to the A8 one-statement-per-line form;
  * equ/log/<name>  parameter file required by ASTRA 8 (A7 had none): copied from a template,
                    with AB/ABC/RTOR/... taken from the first --exp file so ABC <= AB.
--exp FILE copies exp files from <A7>/exp to <A8>/exp unchanged (the patched ASTRA 8 readers
accept the A7 forms used there).
The A8 name of a model is its A7 name with blanks and '+' replaced by '_' and a trailing
'.model' dropped (e.g. "NBI_nubeam_Int_CXRS.model" -> "NBI_nubeam_Int_CXRS").

Existing files in the A8 tree are never overwritten. Formulas present in both trees but with
different text are reported: ASTRA 8 keeps its own version, check them by hand.
"""

import argparse, json, os, re, sys

# A7 NEOCL4 (MAstraN/sbr/neocl4.f) wrote species blocks of 20 into WORK(j,101...300);
# A8 NEOCL4 (sbr/nclass_mod.f90) writes the same quantities into module arrays.
SPECIES = ['e', 'i', 'p', 'd', 't', 'he3', 'he4', 'imp1', 'imp2', 'imp3']
FIELDS = ['gamma_%s_nc', 'dn_%s_nc', 'cn_%s_nc', 'qcond_%s_nc', 'x%s_nc', 'c%s_nc',
          'qen_%s_nc', 'bs_p%s_nc', 'bs_t%s_nc', 'polflow_%s_nc']
WORK_MAP = {101 + 20*js + jf: fld %sp for js, sp in enumerate(SPECIES) for jf, fld in enumerate(FIELDS)}
WORK_MAP.update({301: 'jbs_nc', 302: 'jext_nc', 303: 'cc_nc'})
# 304 (A7: density of main ions) has no exact A8 counterpart (ni_nc sums all thermal ions)

WORK_RE = re.compile(r'\bwork\s*\(\s*j\s*,\s*(\d+)\s*\)', re.IGNORECASE)
INCLUDE_RE = re.compile(r"include\s+'fml/([^']+)'", re.IGNORECASE)
NAME_RE = re.compile(r'\b[A-Za-z][A-Za-z0-9_]*\b')
# Exp scalars carried into equ/log/<model>
LOG_FROM_EXP = ('AB', 'ABC', 'AWALL', 'RTOR', 'ELONG', 'ELONM', 'TRIAN', 'TRICH', 'BTOR', 'IPL',
                'AMJ', 'ZMJ')


def a8_name(a7_name):
    name = a7_name[:-6] if a7_name.endswith('.model') else a7_name
    return re.sub(r'[ +]', '_', name)


def convert_work(text, warnings):
    def repl(m):
        n = int(m.group(1))
        if n in WORK_MAP:
            return '%s(j)' %WORK_MAP[n]
        warnings.append('work(j,%d) has no ASTRA 8 counterpart, left as is' %n)
        return m.group(0)
    return WORK_RE.sub(repl, text)


def fml_a7_to_a8(lines):
    """A7 fixed form (C comments, tab/col-6 continuation) -> A8 free form."""
    out = []
    for raw in lines:
        line = raw.rstrip('\n').rstrip()
        if not line.strip():
            continue
        if line[0] in 'Cc*!':
            out.append('!' + line[1:].replace('\t', ' '))
            continue
        tabbed = line.startswith('\t')
        if tabbed:
            body = line[1:]
            cont = bool(body) and body[0] in '123456789'
            stmt = body[1:] if cont else body
        else:
            cont = len(line) > 5 and line[:5].strip() == '' and line[5] not in ' 0'
            stmt = line[6:] if cont else line
        stmt = stmt.replace('\t', ' ').strip()
        stmt = stmt.split('!')[0].rstrip() if not stmt.lower().startswith('include') else stmt
        if cont and out and not out[-1].startswith('!'):
            out[-1] = out[-1] + ' ' + stmt
        elif stmt:
            out.append(stmt)
    return '\n'.join(out) + '\n'


def include_fnc(text, a8_fml, a8_fnc, warnings, fml):
    """include 'fml/x' of a formula that A8 has only as function fnc/x.f90: drop the include
    and evaluate the function at the current grid point, X -> XR(J*HRO)."""
    for inc in INCLUDE_RE.findall(text):
        name = inc.lower()
        if name in a8_fml or name not in a8_fnc:
            continue
        text = re.sub(r"(?im)^\s*include\s+'fml/%s'\s*\n" %re.escape(inc), '', text)
        body, comments = [], []
        for line in text.splitlines():
            if line.lstrip().startswith('!'):
                body.append(line)
            else:
                body.append(re.sub(r'\b%s\b' %re.escape(name), '%sR(J*HRO)' %name.upper(), line, flags=re.IGNORECASE))
        text = '\n'.join(body) + '\n'
        warnings.append('fml/%s: %s is an A8 function, used as %sR(J*HRO)' %(fml, name.upper(), name.upper()))
    return text


def norm(text):
    return re.sub(r'\s+', '', '\n'.join(l for l in text.lower().splitlines()
                                        if l.strip() and l.lstrip()[0] not in 'c!*'))


def exp_scalars(path):
    vals = {}
    with open(path, errors='replace') as f:
        for line in f:
            line = line.split('!')[0]
            words = line.split()
            if len(words) < 2 or line[0].isspace():
                continue
            key = words[0].upper()
            key = key[:-1] if key.endswith('X') and key[:-1] in LOG_FROM_EXP else key
            if key not in LOG_FROM_EXP or key in vals:
                continue
            # value column 17-22 (old standard "Name Time Value Error"), widened to the word
            j1, j2 = 16, 22
            while j1 > 7 and j1 < len(line) and not line[j1].isspace() and not line[j1-1].isspace():
                j1 -= 1
            while j2 < len(line) and not line[j2-1].isspace() and not line[j2].isspace():
                j2 += 1
            try:
                vals[key] = float(line[j1: j2].strip().upper().replace('D', 'E'))
            except ValueError:
                pass
    return vals


def resolve_a7(path, a7_user):
    """A7 user areas link to the master copy by absolute path (".../ASTRA_7.02/MAstraN/fml/x");
    re-root such links on the directory above <a7 user> when the archive was moved."""
    if os.path.exists(path):
        return path
    if os.path.islink(path):
        target = os.readlink(path)
        root = os.path.dirname(os.path.abspath(a7_user.rstrip('/')))
        parts = target.split('/')
        for k in range(1, len(parts)):
            cand = os.path.join(root, *parts[k:])
            if os.path.exists(cand):
                return cand
    return None


def write_new(path, text, dry, log):
    if os.path.exists(path):
        log.append('kept existing %s' %path)
        return False
    log.append(('would write ' if dry else 'wrote ') + path)
    if not dry:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, 'w') as f:
            f.write(text)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('a7', help='A7 user directory (with equ/, fml/, exp/)')
    ap.add_argument('a8', help='ASTRA 8 installation directory (AWD)')
    ap.add_argument('models', nargs='+', help='model file names in <a7>/equ')
    ap.add_argument('--exp', nargs='*', default=[], help='exp files in <a7>/exp to copy; the first one sets equ/log geometry')
    ap.add_argument('--log-template', default='flux_feqis', help='A8 equ/log file used as template (default flux_feqis)')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    a8_fml = set(os.listdir(os.path.join(args.a8, 'fml')))
    # In A8 many A7 formulas became functions (fnc/<name>.f90, e.g. QOH, WE, ITOT) or built-in
    # variables; an fml with the same name would shadow them, so those are never converted.
    a8_fnc = {f.split('.')[0].lower() for f in os.listdir(os.path.join(args.a8, 'fnc'))}
    with open(os.path.join(args.a8, 'astra_variables.json')) as f:
        a8_vars = {v.lower() for group in json.load(f).values() for v in group}
    builtin = set()
    a7_fml_dir = os.path.join(args.a7, 'fml')
    a7_fml = set(os.listdir(a7_fml_dir)) if os.path.isdir(a7_fml_dir) else set()
    log, warnings, differ = [], [], set()

    geo = exp_scalars(os.path.join(args.a7, 'exp', args.exp[0])) if args.exp else {}
    for exp in args.exp:
        with open(os.path.join(args.a7, 'exp', exp), errors='replace') as f:
            write_new(os.path.join(args.a8, 'exp', exp), f.read(), args.dry_run, log)

    for model in args.models:
        with open(os.path.join(args.a7, 'equ', model), errors='replace') as f:
            text = f.read()
        name = a8_name(model)
        write_new(os.path.join(args.a8, 'equ', name), convert_work(text, warnings), args.dry_run, log)

        # formulas: every name in the model that is an A7 fml, plus their includes
        code = '\n'.join(l.split('!')[0] for l in text.splitlines())   # names in comments do not count
        todo = [w.lower() for w in set(NAME_RE.findall(code)) if w.lower() in a7_fml]
        seen = set()
        while todo:
            fml = todo.pop()
            if fml in seen:
                continue
            seen.add(fml)
            a7_path = resolve_a7(os.path.join(a7_fml_dir, fml), args.a7)
            if a7_path is None:
                if fml not in a8_fml:
                    warnings.append('fml/%s: broken link in the A7 tree, not converted' %fml)
                continue
            with open(a7_path, errors='replace') as f:
                a7_lines = f.readlines()
            todo += [m.lower() for m in INCLUDE_RE.findall(''.join(a7_lines)) if m.lower() in a7_fml]
            new_text = fml_a7_to_a8(a7_lines)
            if re.search(r'\bsystemqq\b', new_text, re.IGNORECASE):
                # Intel-only file operations (e.g. fml/save copying equ/exp into dat/)
                warnings.append('fml/%s calls Intel-only SYSTEMQQ: replaced by a stub %s=0' %(fml, fml.upper()))
                new_text = '! A7 fml/%s used SYSTEMQQ (Intel only); stub written by a7to8.py\n%s=0.\n' %(fml, fml.upper())
            new_text = include_fnc(new_text, a8_fml, a8_fnc, warnings, fml)
            if fml not in a8_fml and (fml in a8_fnc or fml in a8_vars):
                builtin.add(fml)
                continue
            if fml in a8_fml:
                with open(os.path.join(args.a8, 'fml', fml), errors='replace') as f:
                    if norm(f.read()) != norm(new_text):
                        differ.add(fml)
                continue
            write_new(os.path.join(args.a8, 'fml', fml), new_text, args.dry_run, log)

        # equ/log/<model>
        tpl = os.path.join(args.a8, 'equ', 'log', args.log_template)
        with open(tpl) as f:
            tlines = f.readlines()
        g = dict(geo)
        if 'AB' in g:
            g.setdefault('ABC', g['AB'])
            g['ABC'] = min(g['ABC'], g['AB'])
            g.setdefault('AWALL', g['AB'])
        out = []
        for line in tlines:
            key = line.split('=')[0].strip() if '=' in line else ''
            if key in g:
                line = '%-6s = %10.3E\n' %(key, g[key])
            out.append(line)
        write_new(os.path.join(args.a8, 'equ', 'log', name), ''.join(out), args.dry_run, log)

    for line in log:
        print(line)
    for w in sorted(set(warnings)):
        print('WARNING:', w)
    if builtin:
        print('NOTE: built into ASTRA 8 (function or variable), A7 formula not copied:', ' '.join(sorted(builtin)))
    if differ:
        print('NOTE: formulas differ between A7 and A8, A8 version kept:', ' '.join(sorted(differ)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
