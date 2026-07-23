#!/usr/bin/env python3

import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(
        description="Generate vmecinput.dat from template with new CURTOR, PHIEDGE, NEQUIL values."
    )
    parser.add_argument("curtor" , type=float, help="Value to set for CURTOR")
    parser.add_argument("phiedge", type=float, help="Value to set for PHIEDGE")
    parser.add_argument("nequil", type=int, help="Value to set for NEQUIL")

    args = parser.parse_args()

    template_path = Path("dat/vmecinput_template.dat")
    output_path   = Path("dat/vmecinput.dat")

    if not template_path.exists():
        raise FileNotFoundError(f"Template file not found: {template_path}")

    text = template_path.read_text()

    text = text.replace("$CURTOR", str(args.curtor))
    text = text.replace("$PHIEDGE", str(args.phiedge))
    print(args.nequil)
    text = text.replace("$NEQUIL", str(args.nequil))

    output_path.write_text(text)

    print(f"Written updated VMEC input to: {output_path}")


if __name__ == "__main__":
    main()
