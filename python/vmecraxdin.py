#!/usr/bin/env python3
import numpy as np
import os
import argparse
import ast
import re
from vmecmodin import replaceParams, replaceArray


def write_vmecinput(raxis, zaxis, curtor, phiedge, nequil, template_file='vmec_io/vmecinput_template.dat', output_file='dat/vmecinput.dat'):
    """Replace raxis and zaxis in a VMEC template and save to dat/vmecinput.dat"""

# Ensure both arrays are the same length
    if len(raxis) != len(zaxis):
        raise ValueError(f"raxis and zaxis must have the same length, got {len(raxis)} and {len(zaxis)}")

# Read template
    with open(template_file, 'r') as f:
        text = f.read()
# Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Update vmecinpit.dat with new raxis, zaxis values
    text = replaceArray(text, ' RAXIS'   , raxis)
    text = replaceArray(text, ' ZAXIS'   , zaxis)

    print('RAXIS, ZAXIS, curtor, phiedge, nequil', raxis, zaxis, curtor, phiedge, nequil)

# Write modified file
    with open(output_file, 'w') as f:
        f.write(replaceParams(text, curtor, phiedge, nequil))

    print(f"Modified VMEC input file written to {output_file} (length {len(raxis)})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate VMEC input file from template")
    parser.add_argument('raxis', type=str  , help='Radial axis array (e.g., "[1.0, 1.1, 1.2]")')
    parser.add_argument('zaxis', type=str  , help='Vertical axis array (e.g., "[0.0, 0.1, 0.2]")')
    parser.add_argument("curtor"  , type=float, help="Value to set for CURTOR")
    parser.add_argument("phiedge" , type=float, help="Value to set for PHIEDGE")
    parser.add_argument("nequil"  , type=int  , help="Value to set for NEQUIL")
    args = parser.parse_args()

    write_vmecinput(args.raxis, args.zaxis, args.curtor, args.phiedge, args.nequil)
