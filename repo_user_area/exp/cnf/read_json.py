import json
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as patches


def parallPlot(ax, r, z, dr, dz, angle_h, angle):
    angle_h = np.radians(angle_h)
    angle   = np.radians(angle)
    x = [r-0.5*dr-0.5*dz*np.cos(angle), r+0.5*dr-0.5*dz*np.cos(angle), r+0.5*dr+0.5*dz*np.cos(angle), r-0.5*dr+0.5*dz*np.cos(angle)]
    y = [z-0.5*dz, z-0.5*dz, z+0.5*dz, z+0.5*dz]
    ax.add_patch(patches.Polygon(xy=list(zip(x,y)), fill=False))

f_json = 'aug_description_in.json'
#f_json = 'iter_description_in.json'
#f_json = 'jet_description_in.json'

with open(f_json) as fjson:
    in_d = json.load(fjson)

fig = plt.figure(1, (14, 9))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
jbeg = 0
print(len(in_d['Rvessel']))

if 'R_coil' in in_d.keys():
    n_coils = len(in_d['R_coil'])
    for jcoil in range(n_coils):
        parallPlot(ax, in_d['R_coil'][jcoil], in_d['Z_coil'][jcoil], in_d['dR_coil'][jcoil], in_d['dZ_coil'][jcoil], in_d['angh_coil'][jcoil], in_d['ang_coil'][jcoil])

for jlen in in_d['contour_len']:
    ax.plot(in_d['Rvessel'][jbeg: jbeg+jlen], in_d['Zvessel'][jbeg: jbeg+jlen])
    jbeg += jlen

plt.show()
