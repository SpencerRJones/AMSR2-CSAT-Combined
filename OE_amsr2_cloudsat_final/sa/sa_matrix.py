'''

Generates Sa matrix for AMSR2-CSAT combined retrievals.

Sa is of shape [nx x nx] where nx is the number of elements in x,
    the state vector. Contains layers of liquid precip water, ice water,
    cloud water path, first 3 eofs of water vapor, mu, SST, and wind speed.

    Current structure for X:
         1-30: layers of rain water:
        31-60: layers of ice:
           61: LOG10(cloud liquid water)
           62: DSD mu
           63: Snow particle density
        64-66: EOF1-3 of water vapor
           67: SST
           68: Wind Speed


'''

import numpy as np


nx = 66

sa = np.zeros([nx,nx], dtype='f')

#---Uncertainties for rain water
for a in np.arange(0,30):
    sa[a,a] = 2.0

#---Uncertainties for ice water
for a in np.arange(30,60):
    sa[a,a] = 2.0

#---Uncertainty for cloud water path
sa[-6,-6] = 2.

#---Uncertainty for mu:
sa[-5,-5] = 0.5

#---Uncertainty for snow particle density
sa[-4,-4] = 0.04

#---Uncertainties for EOFS
sa[-3,-3] = 1.34     #1.34
#sa[-5,-5] = 2.1
#sa[-4,-4] = 0.00001 #0.91
#sa[-3,-3] = 0.00001 #0.40

#---Uncertainties for SST:
#sa[-2,-2] = 4.
sa[-2,-2] = 0.75**2 ###

#---Uncertainties for wind speed:
sa[-1,-1] = 4.




with open('sa.mtrx', 'wb') as f:
    sa.T.astype('f').tofile(f, sep='', format='unformatted')


