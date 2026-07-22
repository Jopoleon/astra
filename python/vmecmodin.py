#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def replace_parameter(text, name, value):
    """
    Replace the first occurrence of a VMEC-style parameter assignment like:
        NAME = something
    with:
        NAME = <value>
    Preserves leading whitespace and anything after the value (like comments).
    """
    pattern = re.compile(rf"^(\s*{name}\s*=\s*)([^!\n\r]*)(.*)$", re.MULTILINE)

    def repl(match):
        prefix = match.group(1)
        suffix = match.group(3)  # comments or trailing stuff
        return f"{prefix}{value}{suffix}"

    new_text, count = pattern.subn(repl, text, count=1)

    if count == 0:
        raise ValueError(f"Parameter '{name}' not found in template file.")

    return new_text


def main():
    parser = argparse.ArgumentParser(
        description="Generate vmecinput.dat from template with new CURTOR and PHIEDGE values."
    )
    parser.add_argument("curtor", type=float, help="Value to set for CURTOR")
    parser.add_argument("phiedge", type=float, help="Value to set for PHIEDGE")

    args = parser.parse_args()

    template_path = Path("dat/vmecinput_template.dat")
    output_path = Path("dat/vmecinput.dat")

    if not template_path.exists():
        raise FileNotFoundError(f"Template file not found: {template_path}")

    text = template_path.read_text()

    text = replace_parameter(text, "CURTOR", args.curtor)
    text = replace_parameter(text, "PHIEDGE", args.phiedge)

    output_path.write_text(text)

    print(f"Wrote updated VMEC input to: {output_path}")
    print(f"  CURTOR   = {args.curtor}")
    print(f"  PHIEDGE = {args.phiedge}")


if __name__ == "__main__":
    main()
