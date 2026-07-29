import os, re

def parse_fortran_namelist(file_path, group_name):
    params = {}
    in_group = False
    var_re = re.compile(r"^\s*([a-zA-Z0-9_]+)\s*=\s*(?:'([^']*)'|\"([^\"]*)\"|([\d.Ee+-]+))\s*")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Namelist file not found: {file_path}")
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.lower().startswith(f'&{group_name.lower()}'):
                in_group = True
                continue
            if line.startswith('/') and in_group:
                break
            if in_group and not line.startswith('!'):
                m = var_re.match(line)
                if m:
                    key = m.group(1).lower()
                    sq, dq, num = m.group(2), m.group(3), m.group(4)
                    if sq is not None:
                        params[key] = sq
                    elif dq is not None:
                        params[key] = dq
                    elif num is not None:
                        try:
                            params[key] = int(num)
                        except ValueError:
                            params[key] = float(num)
    return params
