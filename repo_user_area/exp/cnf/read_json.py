import json
import numpy as np
import matplotlib.pylab as plt


f_json = 'iter_description_in.json'

with open(f_json) as fjson:
    in_d = json.load(fjson)

plt.figure(1)
plt.subplot(1, 1, 1, aspect='equal')
jbeg = 0
for jlen in in_d['contour_len']:
    plt.plot(in_d['Rvessel'][jbeg: jbeg+jlen], in_d['Zvessel'][jbeg: jbeg+jlen])
    jbeg += jlen
plt.show()
