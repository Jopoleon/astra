#!/usr/bin/env python3
import numpy as np
import os
import argparse
import ast
import re

def write_vmecinput(raxis_cc, zaxis_cs, template_file='vmec_io/vmecinput_template.dat', output_file='dat/vmecinput.dat'):
    """Replace raxis_cc and zaxis_cs in a VMEC template and save to dat/vmecinput.dat"""
    # Convert to numpy arrays
    raxis_cc = np.asarray(raxis_cc, dtype=np.float64)
    zaxis_cs = np.asarray(zaxis_cs, dtype=np.float64)

    # Ensure both arrays are the same length
    if len(raxis_cc) != len(zaxis_cs):
        raise ValueError(f"raxis_cc and zaxis_cs must have the same length, got {len(raxis_cc)} and {len(zaxis_cs)}")

    # Read template
    with open(template_file, 'r') as f:
        template = f.read()

    # Convert arrays to comma-separated strings
    raxis_str = '  '.join(f'{x:.6f}' for x in raxis_cc)
    zaxis_str = '  '.join(f'{x:.6f}' for x in zaxis_cs)

    # Replace placeholders
    filled = re.sub(r'^\s*(RAXIS_CC\s*=).*', r'\1 ' + raxis_str, template, flags=re.MULTILINE)
    filled = re.sub(r'^\s*(ZAXIS_CS\s*=).*', r'\1 ' + zaxis_str, filled, flags=re.MULTILINE)
    filled = re.sub(r'^\s*(RAXIS\s*=).*', r'\1 ' + raxis_str, template, flags=re.MULTILINE)
    filled = re.sub(r'^\s*(ZAXIS\s*=).*', r'\1 ' + zaxis_str, filled, flags=re.MULTILINE)
    filled = filled.replace("$CURTOR" , str(args.curtor))
    filled = filled.replace("$PHIEDGE", str(args.phiedge))
    filled = filled.replace("$NEQUIL" , str(args.nequil))

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Write modified file
    with open(output_file, 'w') as f:
        f.write(filled)

    print(f"Modified VMEC input file written to {output_file} (length {len(raxis_cc)})")

def parse_array(array_str):
    """Parse a string representation of a list into a Python list safely"""
    return ast.literal_eval(array_str)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate VMEC input file from template")
    parser.add_argument('raxis_cc', type=str  , help='Radial axis array (e.g., "[1.0, 1.1, 1.2]")')
    parser.add_argument('zaxis_cs', type=str  , help='Vertical axis array (e.g., "[0.0, 0.1, 0.2]")')
    parser.add_argument("curtor"  , type=float, help="Value to set for CURTOR")
    parser.add_argument("phiedge" , type=float, help="Value to set for PHIEDGE")
    parser.add_argument("nequil"  , type=int  , help="Value to set for NEQUIL")
    args = parser.parse_args()

    raxis_cc = parse_array(args.raxis_cc)
    zaxis_cs = parse_array(args.zaxis_cs)

    write_vmecinput(raxis_cc, zaxis_cs, args.curtor, args.phiedge, args.nequil)
