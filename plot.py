import sys, re
import numpy as np
from matplotlib import pyplot as plt

def remove_all_whitespace(x):
    pattern = re.compile(r'\s+')
    return re.sub(pattern, '', x)

def whitespace_to_semicol(x):
    pattern = re.compile(r'\s+')
    return re.sub(pattern, ';', x)

def split(x):
    return whitespace_to_semicol( x ).split( ';' )[1:-1]

efermi = float( sys.argv[1] )

################################################################################
# parse KPOINTS file
################################################################################

with open('KPOINTS', 'r') as f:
    data = f.readlines()

# assume first line looks like:
# some text: X1-X2-X3-...-XN|Y1-...-YM|Z1-...

path = remove_all_whitespace( data[0].split(':')[1] )

labels = re.split('-', path)
for i, lab in enumerate(labels):
    if lab=='G': labels[i] = '$\Gamma$'

npl = int( remove_all_whitespace( data[1] ) )

################################################################################
# parse EIGENVAL file
################################################################################

with open('EIGENVAL', 'r') as f:
    data = f.readlines()

natoms, _, _, _ = ( int(x) for x in split( data[0] ) )
nocc, nkp, nbands = ( int(x) for x in split( data[5] ) )
pos = 7
bands = np.zeros( (nkp, nbands) )
kpoints = np.zeros( (nkp, 3) )

for ik in range(nkp):
    kpoints[ik] = np.array( [ float(x) for x in split( data[pos] )[:-1] ] )
    pos += 1
    for ib in range(nbands):
        bands[ik,ib] = float( split( data[pos+ib] )[1] ) - efermi
    pos += nbands+1

################################################################################
# plotting
################################################################################

nsegments = int(nkp/npl)

k = np.zeros(nkp)
labelsk = np.zeros( len(labels) )
offset = np.linalg.norm( kpoints[0] )

for i in range(nsegments):
    labelsk[i] = offset
    for j in range(npl):
        k[i*npl+j] = np.linalg.norm( kpoints[i*npl+j] - kpoints[i*npl]) + offset
    #if ( i+1 < nsegments ): 
    offset += np.linalg.norm( kpoints[(i+1)*npl-1] - kpoints[i*npl-1])
labelsk[-1] = offset


for i in range(520, 533):
    if i== 528: color='crimson'
    else: color='royalblue'
    plt.plot(k, bands[:,i], color=color, linewidth=1)

for x in labelsk:
    plt.axvline(x, color='lightgrey', linewidth=2)

plt.hlines(0, min(k), max(k), color='darkorange')
plt.xlim(min(k), max(k))
plt.xticks( labelsk, labels=labels, size=20)
plt.yticks(fontsize=15)
plt.ylabel('$E - E_{F}$ (eV)', size=20)
plt.tight_layout()
plt.show()
