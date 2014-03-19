from math import sqrt
from random import random, randint, gauss

seconds = 1           # if 1, toas are in seconds. 0 = MJD
N = 600               # Number of points
pfract = 0.6          # Number of photons from PSR / Number from background
freq = 34.1234565     # Pulsar Frequency (hz)
#fdot = -2.34556e-10   # Pulsar F-dot (hz/s)
fdot = 0.0            # Pulsar F-dot (hz/s)
T = 20000.0           # Total integration time (s)
width = 0.1           # +/-1 sigma width of Gaussian pulse in phase

numrot = int(T * freq) # approximately
toas = []
if (seconds):
    factor = 1.0
    offset = 0.0
else:
    factor = 86400.0
    offset = 0.0
for point in range(N):
    if (random() < pfract):
        toas.append(randint(0, numrot) + gauss(0.0, width/2.0))
    else:
        toas.append(randint(0, numrot) + random())
toas.sort()
file = open("test.events", "w")
if (fdot == 0.0):
    for point in range(N):
        toas[point] = toas[point] / (freq * factor) + offset
        file.write("%19.10f\n" % toas[point])
else:
    fdot = fdot * 2.0
    for point in range(N):
        toas[point] = ((sqrt(freq * freq + 2 * toas[point] \
                             * fdot) - freq) / fdot) \
                             / factor + offset
        file.write("%19.10f\n" % toas[point])
file.close()
