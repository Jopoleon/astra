#!/usr/bin/env python3

import re

input_file  = "vmec_io/vmecinput.dat"
output_file = "vmec_io/nsurfaces.dat"

ns_array = []

with open(input_file, "r") as f:

    for line in f:

        # Remove comments
        line = line.split("!")[0]

        # Find ns_array line
        if "ns_array" in line.lower():

            # Extract everything after =
            rhs = line.split("=")[1]

            # Extract all integers
            numbers = re.findall(r'\d+', rhs)

            ns_array.extend(int(n) for n in numbers)


if not ns_array:
    raise ValueError("ns_array not found in vmecinput.dat")


# Last value is the actual number of surfaces
nsurfaces = ns_array[-1]


# Write output
with open(output_file, "w") as f:
    f.write(str(nsurfaces) + "\n")


print("ns_array =", ns_array)
print("Number of surfaces =", nsurfaces)
print("Written to", output_file)
