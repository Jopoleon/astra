#!/usr/bin/env python3

import argparse, ast, re
from pathlib import Path
import numpy as np


def replaceArray(text, arrName, arrStr):

    arrList = ast.literal_eval(arrStr)
    arr_str = '  '.join(f'{x:.6f}' for x in arrList)
    text = re.sub(r'^\s*(' + arrName + '\s*=).*', r'\1 ' + arr_str, text, flags=re.MULTILINE)
    return text

    
def replaceParams(text, curtor, phiedge, nequil):

    text = text.replace("$CURTOR" , '%12.4e' %curtor)
    text = text.replace("$PHIEDGE", '%12.4e' %phiedge)
    text = text.replace("$NEQUIL" , '%d'     %nequil)

    return text


def main():

    parser = argparse.ArgumentParser(
        description="Generate vmecinput.dat from template with new CURTOR, PHIEDGE, NEQUIL values."
    )
    parser.add_argument("curtor" , type=float, help="Value to set for CURTOR")
    parser.add_argument("phiedge", type=float, help="Value to set for PHIEDGE")
    parser.add_argument("nequil" , type=int  , help="Value to set for NEQUIL")

    args = parser.parse_args()

    template_path = Path("vmec_io/vmecinput_template.dat")
    output_path   = Path("dat/vmecinput.dat")

    if not template_path.exists():
        raise FileNotFoundError(f"Template file not found: {template_path}")

    text = replaceParams(template_path.read_text(), args.curtor, args.phiedge, args.nequil)

    output_path.write_text(text)

    print(f"Written updated VMEC input to: {output_path}")


if __name__ == "__main__":
    main()
