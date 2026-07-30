#!/usr/bin/env python3

import argparse, re, os
import numpy as np
from parse_fortran_nml import parse_fortran_namelist


def replaceArray(text, arrName, arrStr):
    return re.sub(r'^\s*( ' + arrName + r'\s*=).*', r'\1 ' + arrStr, text, flags=re.MULTILINE)

    
def replaceParams(text, curtor, phiedge, nequil):

    text = text.replace("$CURTOR" , '%12.4e' %curtor)
    text = text.replace("$PHIEDGE", '%12.4e' %phiedge)
    text = text.replace("$NEQUIL" , '%d'     %nequil)

    return text


def write_vmecinput(curtor, phiedge, nequil, raxis=None, zaxis=None, template_file='vmec_io/vmecinput_template.dat', namelist_path='vmec_io/stell_files.nml'):
    """Replace CURTOR, PHIEDGE, NEQUIL (optional RAXIS, ZAXIS) from vmec_io/vmecinput_template.dat into dat/vmecinput.dat"""

# Read template
    if not os.path.isfile(template_file):
        raise FileNotFoundError(f"Template file not found: {template_path}")

    nl_vmec = parse_fortran_namelist(namelist_path, 'VMEC_TO_ASTRA_INPUTS')
    print(nl_vmec.keys())
    output_file = nl_vmec['vmec_in_file']

    with open(template_file, 'r') as f:
        text = f.read()

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    if raxis is not None:
        raxis = raxis.strip()
        zaxis = zaxis.strip()
        if len(raxis.split()) != len(zaxis.split()):
            raise ValueError(f"raxis and zaxis must have the same length, got {len(raxis)} and {len(zaxis)}")
        text = replaceArray(text, 'RAXIS', raxis)
        text = replaceArray(text, 'ZAXIS', zaxis)

# Write modified file, replacing CURTOR, PHIEDGE, NEQUIL
    with open(output_file, 'w') as f:
        f.write(replaceParams(text, curtor, phiedge, nequil))

    print(f"Written updated VMEC input to: {output_file}")


def main():

    parser = argparse.ArgumentParser(
        description="Generate vmecinput.dat from template with new CURTOR, PHIEDGE, NEQUIL values."
    )
    parser.add_argument("curtor" , type=float, help="Value to set for CURTOR")
    parser.add_argument("phiedge", type=float, help="Value to set for PHIEDGE")
    parser.add_argument("nequil" , type=int  , help="Value to set for NEQUIL")
    parser.add_argument('-r', '--raxis', type=str, help='Radial axis array (e.g., "[1.0, 1.1, 1.2]")'  , required=False)
    parser.add_argument('-z', '--zaxis', type=str, help='Vertical axis array (e.g., "[0.0, 0.1, 0.2]")', required=False)

    args = parser.parse_args()

    write_vmecinput(args.curtor, args.phiedge, args.nequil, raxis=args.raxis, zaxis=args.zaxis)


if __name__ == "__main__":
    main()
