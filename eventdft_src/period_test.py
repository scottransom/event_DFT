from TableIO import readColumns
from presto import spectralpower
from periodogram import *
from time import clock

file = 'rosat_crab_3k.toa'
toas = ravel(readColumns(file, '#'))
toas = (toas - toas[0]) * 86400.0
T = toas[-1] - toas[0]
df = 1.0 / (10.0 * T)

p0 = 0.033468239225/2
numf = 100000
lof = 1.0 / p0 - df * numf / 2.0
freqs = arange(numf) * df + lof
beg = clock()
mypows = spectralpower(toafft(toas, len(toas), lof, df, numf))
end = clock()
print '     TOAFFT took ', end-beg, ' sec.  (%d freqs/sec)' % (numf/(end-beg))
beg = clock()
pows = periodogram(ones(len(toas), 'd'), toas, len(toas), lof, df, numf)
end = clock()
print 'periodogram took ', end-beg, ' sec.  (%d freqs/sec)' % (numf/(end-beg))

def toafft(toas, freqs, days=0):
    if days:
        toas = (toas - toas[0]) * 86400.0
    else:
        toas = toas - toas[0]
    twid = zeros(len(freqs), 'D')
    twid.imag = -2.0 * pi * freqs
    return spectralpower(add.reduce(exp(\
        multiply.outer(toas, twid)))) / len(toas)

