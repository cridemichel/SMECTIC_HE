#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator)
#plt.title('Dati')
ax = plt.gca()
#ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.02))
plt.xlabel(r'$\phi$', size=16)
plt.ylabel(r'$\beta P v_0$', size=16)
x, y = np.loadtxt('eos_iso.dat', comments=['#'], usecols=(0,1), unpack=True)
plt.plot(x, y, 'go--', label='isotropic')
x, y = np.loadtxt('eos_lamellar.dat', comments=['#'], usecols=(0,1), unpack=True)
plt.plot(x, y, 'ro--', label='lamellar')
x, y = np.loadtxt('eos_glass.dat', comments=['#'], usecols=(0,1), unpack=True)
plt.plot(x, y, 'co--', label='glass')
plt.legend()
plt.savefig('eos_polydisp.png')
plt.show()
