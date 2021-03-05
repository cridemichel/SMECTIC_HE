#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter)
argv =sys.argv
if len(argv) < 2:
    print("Please supply a file with data to plot!")
    quit()
fn=argv[1]
#plt.title('Dati')
ax = plt.gca()
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.02))
plt.xlabel(r'$\phi$', size=16)
plt.ylabel(r'$\beta P v_0$', size=16)
x, y = np.loadtxt(fn, comments=['#'], usecols=(0,1), unpack=True)
plt.plot(x, y, 'ro--')
plt.savefig('eos_polydisp.png')
plt.show()
