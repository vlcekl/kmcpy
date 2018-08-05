#!//anaconda/envs/py36/bin/python
#
# File name:   mklatt.py
# Date:        2018/08/02 16:15
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np

def make_sc(box):
    print('sc', *box)
    return

def make_bcc(box):
    print('bcc', *box)
    return

def make_fcc(box):

    # box dimensions in lattice units
    lx, ly, lz = box
    print(lx*ly*lz)
    print('fcc', 2*lx, ly, lz)

    # layer number
    for iz in range(lz):
        # layer structure
        for iy in range(ly):
            for ix in range(lx):
                rx = 2*ix + (iy + iz)%2
                print('Ni', rx, iy, iz)

    return

if __name__ == "__main__":

    # dictionary of lattice build functions
    make_lattice = {'fcc':make_fcc, 'bcc':make_bcc, 'sc':make_sc}

    # read lattice type and parameters
    with open(sys.argv[1], 'r') as f:

        # lattice type
        ltype = re.findall('\S+', f.readline())[0]

        # box dimensions
        box = [int(d) for d in re.findall('\S+', f.readline())]

    # make lattice
    make_lattice[ltype](box)

# end of mklatt.py 
