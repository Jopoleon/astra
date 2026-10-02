#!/usr/bin/env python3
"""result_to_astra for the STRAHL of the ASTRA 7 package (written for ASTRA 8 sbr/a2strahl.f90).

ASTRA 8 runs "result_to_astra <first element>" in the STRAHL instance directory and reads
stdout (results.txt): Nr, then labelled blocks (label = one token, then Nr values):
rho_pol, Zeff, Prad_tot [W/m^3], n_main [m^-3], Prad_<sp> [W/m^3] for each species,
Prad_main [W/m^3], n_<sp> [m^-3], Zavg_<sp>, Dneo_<sp> [m^2/s], Vneo_<sp> [m/s],
nesrc_<sp> [m^-3/s].

The old STRAHL writes one netCDF file per element: result/<El>strahl_result.dat (El = 2 chars, e.g. W_).
Used here (last time slice): rho_poloidal_grid, electron_density, impurity_density
(per charge state), impurity_radiation (index nion+3 = total, index nion = main-ion
bremsstrahlung, see STRAHL emissiv). Main-ion density and Zeff follow from
quasi-neutrality with all species. Not available in this STRAHL output: neoclassical
D/V averages and the electron source from impurity ionisation -> written as zeros.
"""
import re, sys, glob, os
import numpy as np
import netCDF4

def species_and_main():
    txt = open('param_files/sparams.dat').read().split('\n')
    zmain, elems = 1.0, []
    for i, ln in enumerate(txt):
        if re.match(r'cv +background ion', ln):
            zmain = float(txt[i + 1].split()[1])
        if re.match(r'cv +number of impurities', ln):
            n = int(txt[i + 1].split()[0])
            j = i + 2
            while not re.match(r'cv +element', txt[j]):
                j += 1
            elems = [txt[j + 1 + k].split()[0] for k in range(n)]
            break
    return zmain, elems

def main():
    zmain, elems = species_and_main()
    if not elems and len(sys.argv) > 1:
        elems = [sys.argv[1]]
    data = []
    for el in elems:
        cand = [os.path.join('result', el + 'strahl_result.dat'),
                os.path.join('result', el + '_strahl_result.dat')]
        d = netCDF4.Dataset([c for c in cand if os.path.exists(c)][0])
        nz = np.array(d.variables['impurity_density'][-1], dtype=float) * 1e6     # (nion, nr) m^-3
        rad = np.array(d.variables['impurity_radiation'][-1], dtype=float) * 1e6  # (nion+3, nr) W/m^3
        nion = nz.shape[0]
        rec = dict(
            rpol=np.array(d.variables['rho_poloidal_grid'][:], dtype=float),
            ne=np.array(d.variables['electron_density'][-1], dtype=float) * 1e6,
            nz=nz, z=np.arange(nion, dtype=float),
            prad_imp=rad[nion + 2] - rad[nion - 1],
            prad_main=rad[nion - 1],
        )
        d.close()
        data.append(rec)
    rpol, ne = data[0]['rpol'], data[0]['ne']
    nr = len(rpol)
    sum_zn = sum((r['z'][:, None] * r['nz']).sum(0) for r in data)
    sum_z2n = sum((r['z'][:, None] ** 2 * r['nz']).sum(0) for r in data)
    nmain = np.maximum(ne - sum_zn, 0.0) / zmain
    zeff = (zmain ** 2 * nmain + sum_z2n) / np.maximum(ne, 1e-30)
    pradmain = data[0]['prad_main']
    pradtot = pradmain + sum(r['prad_imp'] for r in data)

    def block(label, arr):
        print(label)
        print(' '.join('%.6e' % (v if abs(v) > 1e-30 else 0.0) for v in arr))

    print(nr)
    block('rho_pol', rpol)
    block('Zeff', zeff)
    block('Prad_tot', pradtot)
    block('n_main', nmain)
    for el, r in zip(elems, data):
        block('Prad_' + el, r['prad_imp'])
    block('Prad_main', pradmain)
    for el, r in zip(elems, data):
        block('n_' + el, r['nz'].sum(0))
    for el, r in zip(elems, data):
        ntot = r['nz'].sum(0)
        block('Zavg_' + el, np.where(ntot > 0, (r['z'][:, None] * r['nz']).sum(0) / np.maximum(ntot, 1e-30), 0.0))
    for el in elems:
        block('Dneo_' + el, np.zeros(nr))
    for el in elems:
        block('Vneo_' + el, np.zeros(nr))
    for el in elems:
        block('nesrc_' + el, np.zeros(nr))

if __name__ == '__main__':
    main()
