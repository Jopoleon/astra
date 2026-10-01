import json
import numpy as np
import matplotlib.pylab as plt


def parallPlot(r, z, dr, dz, angle):

    angle   = np.radians(angle)
    dr_2 = 0.5*dr
    dz_2 = 0.5*dz
    dz_angle = dz_2/np.tan(angle)

    x1 = r - dr_2 - dz_angle
    x2 = r + dr_2 - dz_angle
    x3 = r + dr_2 + dz_angle
    x4 = r - dr_2 + dz_angle
    x = [x1, x2, x3, x4, x1]
    y1 = z - dz_2
    y2 = z + dz_2
    y = [y1, y1, y2, y2, y1]
    plt.plot(x, y, color='darkgrey')


if __name__ == '__main__':
    
    f_json = 'aug_description_in.json'
#    f_json = 'iter_description_in.json'
#    f_json = 'jet_description_in.json'

    with open(f_json) as fjson:
        in_d = json.load(fjson)

    plt.figure(1, (9, 10))
    plt.subplot(1, 1, 1, aspect='equal')
    jbeg = 0

    for jcoil, jlen in enumerate(in_d['contour_len']):
        if in_d['contour_color'][jcoil] > 0:
            plt.plot(in_d['Rvessel'][jbeg: jbeg+jlen], in_d['Zvessel'][jbeg: jbeg+jlen], 'b-')
        jbeg += jlen

    if 'R_coil' in in_d.keys():
        n_coils = len(in_d['R_coil'])
        for jcoil in range(n_coils):
            plt.plot(in_d['R_coil'][jcoil], in_d['Z_coil'][jcoil], 'go')
            parallPlot(in_d['R_coil'][jcoil], in_d['Z_coil'][jcoil], in_d['dR_coil'][jcoil], in_d['dZ_coil'][jcoil], in_d['ang_coil'][jcoil])

    plt.xlabel('R [m]')
    plt.ylabel('z [m]')
    plt.show()
