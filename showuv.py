#!/usr/bin/python

import numpy as np
from matplotlib import *
import matplotlib.pyplot as plt
import sys

if len(sys.argv) < 2:
    print './showuv.py <output_file>'
    sys.exit(0)

arr = np.loadtxt(sys.argv[1], dtype = float)

rc('font', size = 17)
width, height = rcParams['figure.figsize']
sz = min(width, height)
plt.figure(figsize = (sz, sz))

plt.xlabel('u')
plt.ylabel('v')

u = arr[:, 0]
v = arr[:, 1]

plt.scatter( u,  v, s = 1, c = 'k', edgecolors = 'none')
plt.scatter(-u, -v, s = 1, c = 'r', edgecolors = 'none')

maxu = np.amax(np.abs(u))
maxv = np.amax(np.abs(v))
uvmax = max(maxu, maxv) * 1.05

plt.ylim(-uvmax, uvmax)
plt.xlim(uvmax, -uvmax)

plt.savefig('out.png')
