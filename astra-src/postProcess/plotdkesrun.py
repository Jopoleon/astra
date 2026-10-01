import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def parse_dkes_file(filename):
    """
    Parses DKES database file and returns DataFrame with:

    surface cmul efield L11 L13 L33
    r R B iota xkn ftrap Bsq
    eps_eff g11_ft efield_u g11_er ex_er
    """

    data = []

    current_surface = 0
    geom = None
    cfit = None

    with open(filename, "r") as f:

        for line in f:
            line = line.strip()

            if not line or line.startswith("cc"):
                continue

            # ---------------------------------
            # Geometry line
            # ---------------------------------
            if "r,R,B,io,xkn,ft,<b^2>" in line:
                parts = line.split()

                current_surface += 1

                geom = {
                    "r": float(parts[0]),
                    "R": float(parts[1]),
                    "B": float(parts[2]),
                    "iota": float(parts[3]),
                    "xkn": float(parts[4]),
                    "ftrap": float(parts[5]),
                    "Bsq": float(parts[6]),
                }

                continue

            # ---------------------------------
            # CFIT line
            # ---------------------------------
            if line.startswith("cfit"):
                parts = line.split()

                cfit = {
                    "eps_eff": float(parts[1]),
                    "g11_ft": float(parts[2]),
                    "efield_u": float(parts[3]),
                    "g11_er": float(parts[4]),
                    "ex_er": float(parts[5]),
                }

                continue

            # ignore comment lines
            if line.startswith("c") or line.startswith("e") or line.startswith(">3"):
                continue

            # ---------------------------------
            # Transport rows
            # ---------------------------------
            parts = line.split()

            if len(parts) >= 5:
                try:

                    cmul = np.log(float(parts[0]))
                    efield = float(parts[1])
                    L11 = float(parts[2])
                    L13 = float(parts[3])
                    L33 = float(parts[4])

                    row = {
                        "surface": current_surface,
                        "cmul": cmul,
                        "efield": efield,
                        "L11": L11,
                        "L13": L13,
                        "L33": L33
                    }

                    if geom:
                        row.update(geom)

                    if cfit:
                        row.update(cfit)

                    data.append(row)

                except ValueError:
                    continue

    df = pd.DataFrame(data)

    return df

def save_output_file(df, filename="output.txt"):
    """
    Saves the DataFrame to a text file in the format:
    surface cmul efield L11 L13 L33
    """
    df.to_csv(filename, sep=' ', index=False, header=False, float_format="%.6E")
    print(f"Data saved to {filename}")

def plot_l_elements(df):
    """
    Plots L11, L13, L33 vs cmul for different efield values, separated by surface
    """
    for surface in df['surface'].unique():
        df_surf = df[df['surface'] == surface]
        for L_element in ["L11", "L13", "L33"]:
            fig, ax = plt.subplots()
            for efield in df_surf['efield'].unique():
                df_ef = df_surf[df_surf['efield'] == efield]
                ax.plot(df_ef['cmul'], df_ef[L_element], marker='o', label=f"Er={efield}")
            ax.set_xlabel("cmul")
            ax.set_ylabel(L_element)
            ax.set_title(f"Surface {surface} - {L_element} vs cmul")
            ax.legend()
            plt.savefig(f"{L_element}_surface{surface}.png")
            plt.close(fig)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_dkes.py <dkes_output_file>")
        sys.exit(1)

    filename = sys.argv[1]
    df = parse_dkes_file(filename)
    if df.empty:
        print("No data parsed. Check file format.")
        sys.exit(1)

    # Save output.txt
    save_output_file(df)

    # Make plots
    plot_l_elements(df)
    print("Plots saved for each surface (L11, L13, L33).")
