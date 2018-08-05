#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import re
import numpy as np
from itertools import product

class KMCMove:
    """Class managing kmc moves and event modifications"""

    def __init__(self, latt_type):
        self.latt_type = latt_type
        self.__setup_neighbors()

    def __setup_neighbors(self):
        """Create lists of neighboring (nn & nnn) sites"""

        if self.latt_type == 'fcc':
            nbrlist = []
            # NNN
            nbrlist.append(np.array([ 2, 0, 0]))
            nbrlist.append(np.array([-2, 0, 0]))
            nbrlist.append(np.array([ 0, 2, 0]))
            nbrlist.append(np.array([ 0,-2, 0]))
            nbrlist.append(np.array([ 0, 0, 2]))
            nbrlist.append(np.array([ 0, 0,-2]))
            # NN
            nbrlist.append(np.array([ 1, 1, 0]))
            nbrlist.append(np.array([ 1,-1, 0]))
            nbrlist.append(np.array([-1, 1, 0]))
            nbrlist.append(np.array([-1,-1, 0]))
            nbrlist.append(np.array([ 1, 0, 1]))
            nbrlist.append(np.array([ 1, 0,-1]))
            nbrlist.append(np.array([-1, 0, 1]))
            nbrlist.append(np.array([-1, 0,-1]))
            nbrlist.append(np.array([ 0, 1, 1]))
            nbrlist.append(np.array([ 0, 1,-1]))
            nbrlist.append(np.array([ 0,-1, 1]))
            nbrlist.append(np.array([ 0,-1,-1]))

        self.nbrlist = nbrlist

    def make_lattice(self, xyz, box):
    
        # create simulation box (fill with -1, e.g., no atoms)
        latt = -np.ones((tuple(box)), dtype=int)
    
        # fill box with atoms
        for i, r in enumerate(xyz):
            latt[tuple(r)] = i
    
        # set -2 to empty sites available for adsorption
        # cycle through empty sites and find those neighboring surface atoms
        for r in product(range(box[0]), range(box[1]), range(box[2])):
            if latt[r] != -1:
                continue
    
            # cycle over nearest neighbors
            neighbors = 0
            for dr in self.nbrlist: 
                rj = (r + dr) % box
                if latt[tuple(rj)] > -1:
                    neighbors += 1
    
            if neighbors > 2:
                latt[r] = -2
    
        self.latt = latt
        self.xyz = xyz
        self.box = box

    def make_event_list(self, rates):
 
        event_list = []

        # deposition events
        for r in product(range(self.box[0]), range(self.box[1]), range(self.box[2])):
            if self.latt[r] == -2:
                event = {}
                event['atom'] = None
                event['current'] = None
                event['new'] = r
                event['type'] = 0
                event['rate'] = rates[0]
                event_list.append(event)

        # diffusion events
        for i, ri in enumerate(self.xyz):
            # cycle over nn and nnn sites
            for dr in self.nbrlist:
                rj = (ri + dr) % self.box
                if self.latt[r] == -2:
                    event = {}
                    event['atom'] = i
                    event['current'] = ri
                    event['new'] = rj
                    event['type'] = 1
                    event['rate'] = rates[1]
                    event_list.append(event)
 
        self.event_list =  event_list
        print('event_list', len(event_list))

        return event_list

    def move(self, event):
        """Perform kmc move given by an event (depsotions or diffusion)"""

        if event['type'] == 0: # deposition
            # create a new atom on a lattice
            self.xyz.append(event['new'])
            latt[tuple(self.xyz[-1][:3])] = len(self.xyz)-1 # atom number
        #else: # diffusion
            
        return new_state

        

class EventTree:
    """Class maintaining a binary tree for random event lookup"""

    def __init__(self):
        self.t = []

    def build_tree(self, events):
        # save events internally
        self.events = events
        print('events', len(events))

        # create an array of rates - level 0 - bottom
        psum = [e['rate'] for e in events]
        if len(psum) %2 == 1: # make sure the list has even number of elements
            psum.extend([0.0])
        self.t.append(np.array(psum))

        # create partial summs (iteratively) up to the 2nd highest level
        while len(psum) > 2:
            psum = [psum[i]+psum[i+1] for i in range(0, len(psum), 2)]
            if len(psum) %2 == 1:
                psum.extend([0.0])
            self.t.append(np.array(psum))

        # create top level = sum of all rates
        self.t.append(np.array(sum(psum)))

        print('tree', self.t)

    def get_topnode(self):
        print('topnode', self.t[-1])
        return self.t[-1]

    # print the 'tree' values by levels
    def print_tree(self):
        print('Number of levels (K):', len(self.t))
        print('Number of events:', len(self.t[0]))
        for i, level in enumerate(self.t):
            print('Level:', i)
            print(*level)

    def update_tree(self, new_events):
        """Update tree after an event"""


    def rebalance_tree(self, atoms):
        """
            Remove events related to selected atoms (e.g., layer is filled)
            and rebuild the event tree
        """

        # compile a new list of events
        self.events = [e for e in self.events if e['atom'] not in atoms]

        # rebuild the tree
        self.build_tree(new_events)


    def find_event(self, q):
        """Find and return an event"""

        # cycle through levels (top->down)
        # start with top-level child (k-2) end with level above bottom (1)
        j = 0
        for k in range(len(self.t)-2, 0, -1):
            # left child value
            left = self.t[k][j]

            if q < left:
                j = 2*j
            else:
                q -= left
                j = 2*j + 1
        
        # bottom level - return selected event
        if q < self.t[0][j]:
            return self.events[j]
        else:
            return self.events[j+1]


def read_cfg(file_name):
    """Read input configuration and create simulation box"""

    with open(file_name, 'r') as f:

        # number of atoms
        nat = int(re.findall('\S+', f.readline())[0])

        sarr = re.findall('\S+', f.readline())

        # lattice type
        lat_type = sarr[0]

        # simulation box size
        box = np.array([float(l) for l in sarr[1:]], dtype=int)

        # read particle coordinates and grain identity
        xyz = []
        for i in range(nat):
            xyz.append(np.array([int(l) for l in re.findall('\S+', f.readline())[1:]]))

    return xyz, box

def read_pars(file_name):
    """Read rates for different event types"""

    with open(file_name, 'r') as f:
        rates = []
        for line in iter(f.readline, ''):
            rates.append(float(re.findall('\S+', line)[1]))

    return np.array(rates)


if __name__ == "__main__":

    # read simulation parameters
    with open(sys.argv[1], 'r') as f:

        # time limit (ms)
        t_max = float(re.findall('\S+', f.readline())[1])

        # input configuration file
        incfg_file = re.findall('\S+', f.readline())[1]

        # kmx parameter file - rates (1/ms)
        pars_file = re.findall('\S+', f.readline())[1]


    # read input configuration file and make lattice
    xyz, box = read_cfg(incfg_file)
    kmc = KMCMove('fcc')
    kmc.make_lattice(xyz, box)

    # read read kmx parameters
    rates = read_pars(pars_file)

    # make event list (e.g., identify deposition sites)
    event_list = kmc.make_event_list(rates)

    # Create an initial event tree
    etree = EventTree()
    etree.build_tree(event_list)

    # start main loop
    np.random.seed(42)
    t = 0.0
    while t < t_max:
        # sum of event rates
        Rs = etree.get_topnode()
        print(Rs)

        # generate a random number [0,Rs)
        q = Rs*np.random.random()

        # find next event to occurr
        event_i = etree.find_event(q)
        print('event', event_i)

        # perform event and get new events
        new_state = move(event_i)

        # update binary search tree
        etree.update_tree(new_state)

        # if deposition event, check if a new layer is complete
        if event_i['type'] == 0:
            atoms = check_layer()
            if len(atoms) > 0:
                etree.rebalance_tree(atoms)

        # advance time
        t += -np.log(np.random.random())/Rs

# end of kmc_pld.py 
