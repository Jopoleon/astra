import re


def parse_fortran_value(value):
    """Convert a Fortran namelist value to a Python value."""
    value = value.strip()

    # Character string
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
        return value[1:-1]

    # Logical
    if value.lower() == ".true.":
        return True
    if value.lower() == ".false.":
        return False

    # Fortran double-precision exponent
    value = value.replace("D", "E").replace("d", "e")

    try:
        return int(value)
    except ValueError:
        return float(value)


def parse_fortran_namelist(file_path, group_name):
    """Parse scalar and simple comma-separated values from a Fortran namelist."""
    params = {}
    in_group = False

    var_re = re.compile(
        r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)"
    )

    with open(file_path) as f:
        for line in f:
            line = line.strip()

            # Find requested namelist group
            if not in_group:
                if line.lower().startswith(f"&{group_name.lower()}"):
                    in_group = True
                continue

            # End of namelist group
            if line.startswith("/"):
                break

            # Ignore empty lines and comments
            if not line or line.startswith("!"):
                continue

            match = var_re.match(line)
            if not match:
                continue

            key, value = match.groups()

            # Remove inline comment
            value = value.split("!", 1)[0].strip()

            # Split comma-separated values
            values = [v.strip() for v in value.split(",") if v.strip()]

            parsed = [parse_fortran_value(v) for v in values]

            # Scalar if only one value, otherwise list
            params[key.lower()] = (
                parsed[0] if len(parsed) == 1 else parsed
            )

    return params
