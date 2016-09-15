# ADW, first part of plot_prms.py, which has Q_s calculations as well
# Written 12/2015, modified 9/2016

import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime as dt

prmsout = np.genfromtxt('ToroOut.statvar', delimiter=' ', skip_header=36)
prmsout = prmsout[:-20]
Q = prmsout[:,8] * 0.0283168 # m3/s

date = []
for row in prmsout:
  date.append(dt(int(row[1]), int(row[2]), int(row[3])))

plt.plot(date, Q) # although semilogy can be easier to read
plt.xlabel('Date', fontsize=16)
plt.ylabel('Discharge [m$^3$ s$^{-1}$]', fontsize=16)
plt.show()

