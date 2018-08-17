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
from collections import Counter

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
            nbrlist.append(np.array([ 0, 0,-2]))
            # NN
            nbrlist.append(np.array([ 1, 1, 0]))
            nbrlist.append(np.array([ 1,-1, 0]))
            nbrlist.append(np.array([-1, 1, 0]))
            nbrlist.append(np.array([-1,-1, 0]))
            nbrlist.append(np.array([ 1, 0,-1]))
            nbrlist.append(np.array([-1, 0,-1]))
            nbrlist.append(np.array([ 0, 1,-1]))
            nbrlist.append(np.array([ 0,-1,-1]))
            nbrlist.append(np.array([ 1, 0, 1]))
            nbrlist.append(np.array([-1, 0, 1]))
            nbrlist.append(np.array([ 0, 1, 1]))
            nbrlist.append(np.array([ 0,-1, 1]))
            nbrlist.append(np.array([ 0, 0, 2]))
        else:
            raise ValueError(f'Chosen {self.latt_type} lattice. Currently only FCC lattice is supported.')

        self.nbrlist = nbrlist

    def make_lattice(self, xyz, box):
    
        # create simulation box (fill with -10, e.g., not and FCC site)
        latt = -10*np.ones((tuple(box)), dtype=int)

        # cycle through lattice sites and mark FCC sites with -1
        for r in product(range(box[0]), range(box[1]), range(box[2])):
            if sum(r) % 2 == 0:
                latt[tuple(r)] = -1
    
        # fill box with atoms
        for i, r in enumerate(xyz):
            assert sum(r) % 2 == 0, f'{r} Atom not on FCC lattice!'
            latt[tuple(r)] = i
    
        self.latt = latt
        self.xyz = xyz
        self.box = box
        self.grain = [0 for i in range(len(xyz))]

    def make_event_list(self, rates):
 
        event_list = []

        # Deposition event
        event = {}
        event['type'] = 0
        # single rate times number of surface atomic sites
        event['rate'] = rates[0]*self.box[0]*self.box[1]
        event['atom'] = None
        event['current'] = None
        event['new'] = None

        event_list.append(event)

        # diffusion events
        for i, ri in enumerate(self.xyz):
            # cycle over neighbor sites
            for dr in self.nbrlist[:13]:

                rj = (ri + dr) % self.box

                if self.latt[tuple(rj)] == -1:
                    event = {}
                    event['atom'] = i
                    event['current'] = ri
                    event['new'] = rj
                    event['type'] = 1
                    event['rate'] = rates[1]
                    event_list.append(event)

        self.event_list =  event_list
        #print('event_list', len(event_list))

        return event_list

    def move(self, event, i):
        """
        Perform kmc move given by an event (depsotions or diffusion)
        and update event list
        """

        if event['type'] == 0: # deposition
            search = True
            while search:
                # choose x-y position
                ix = np.random.randint(self.box[0])
                iy = np.random.randint(self.box[1])
                # find iz position
                for iz in range(self.box[2]):
                    if self.latt[ix, iy, iz] == -1:
                        ri = np.array([ix, iy, iz], dtype=int)
                        break
                else:
                    continue

                # search neighbors (needs at least 3)
                neighbors = 0
                grain_numbers = []
                for dr in self.nbrlist:
                    rj = (ri + dr) % self.box
                    iatom = self.latt[tuple(rj)]
                    if iatom > -1:
                        neighbors += 1
                        grain_numbers.append(self.grain[iatom])

                if neighbors > 2:
                    # end search
                    search = False

            # create a new atom
            self.xyz.append(np.array((ix, iy, iz))) 
            # put it on a lattice
            self.latt[tuple(self.xyz[-1])] = len(self.xyz)-1 # atom number to lattice

            # find and adopt the biggest neighboring grain. If none, assign new
            grain_counts = Counter(grain_numbers)
            if 0 in grain_counts:
                del grain_counts[0]

            if any(grain_counts):
                g_number = sorted(grain_counts, key=grain_counts.__getitem__)[-1]
            else:
                g_number = max(self.grain) + 1

            # grain number for each atom
            self.grain.append(g_number)

            print('top atom', ri, ix, iy, iz)

        else: # diffusion
            # cycle over list of surface atoms and their diffusion events
            # choose x-y position
            ix = np.random.randint(self.box[0])
            iy = np.random.randint(self.box[1])
            # find iz position
            for iz in range(self.box[2]-1, 0-1, -1):
                iatom = self.latt[ix, iy, iz]
                if iatom > -1:
                    ri = np.array([ix, iy, iz], dtype=int)
                    break

            # search neighbors (needs at least 3)
            neighbors = 0
            diffusion_directions = []
            for dr in self.nbrlist[:13]:
                rj = (ri + dr) % self.box
                iatom = self.latt[tuple(rj)]
                if iatom == -1:
                    neighbors += 1
                    diffusion_directions.append(dr)

            if len(diffusion_directions) > 0:
                print('diff', diffusion_directions)
                #diff_dir = np.random.choice(np.array(diffusion_directions)) 

                diff_dir = diffusion_directions[np.random.randint(len(diffusion_directions))]
                self.latt[ix, iy, iz] = -1
                new_ri = (ri + diff_dir) % self.box
                self.xyz[iatom] = new_ri
                self.latt[tuple(new_ri)] = iatom

                # find and adopt the biggest neighboring grain. If none, assign new
                grain_numbers = []
                for dr in self.nbrlist[:13]:
                    rj = (new_ri + dr) % self.box
                    jatom = self.latt[tuple(rj)]
                    if jatom > -1:
                        grain_numbers.append(self.grain[iatom])
    
                grain_counts = Counter(grain_numbers)
                if 0 in grain_counts:
                    del grain_counts[0]
                
                if any(grain_counts):
                    g_number = sorted(grain_counts, key=grain_counts.__getitem__)[-1]
                else:
                    g_number = max(self.grain) + 1
                
                # grain number for each atom
                self.grain[iatom] = g_number

        # for atom i in final position ri:
        # change deposition site event for ri
        # find neighboring atoms or deposition sites
        #for dr in self.nbrlist:
        #    rj = (ri + dr) % self.box
        #    if self.latt[r] == -2:
        #        event = {}
        #        event['atom'] = i
        #        event['current'] = ri
        #        event['new'] = rj
        #        event['type'] = 1
        #        event['rate'] = rates[1]
        #        event_list.append(event)

        # new atom position effect on events

        old_events = []
        new_events = []

        return old_events, new_events

    def get_conf(self):
        return self.xyz, self.box, self.grain


class EventTree:
    """Class maintaining a binary tree for random event lookup"""

    def __init__(self):
        self.t = []

    def build_tree(self, events):
        # save events internally
        self.events = events
        #print('events', len(events))

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

        #print('tree', self.t)

    def get_topnode(self):
        #print('topnode', self.t[-1])
        return self.t[-1]

    # print the 'tree' values by levels
    def print_tree(self):
        print('Number of levels (K):', len(self.t))
        print('Number of events:', len(self.t[0]))
        for i, level in enumerate(self.t):
            print('Level:', i)
            print(*level)

    def update_tree(self, new_events, old_events):
        """Update tree after an event"""
        pass


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
            return self.events[j], j
        else:
            return self.events[j+1], j+1


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

def write_cfg(file_name, xyz, box, grain):
    """Write output configuration to xyz file"""

    with open(file_name, 'w') as f:
        # number of atoms and box
        nat = len(xyz)
        f.write(f'{nat}\n')
        f.write(f'fcc {box[0]} {box[1]} {box[2]}\n')

        # write particle coordinates and grain identity
        for r, g in zip(xyz, grain):
            f.write(f'{g} {r[0]} {r[1]} {r[2]}\n')

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
    kmc = KMCMove('bcc')
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
        #print(Rs)

        # generate a random number [0,Rs)
        q = Rs*np.random.random()

        # find next event to occurr
        event, i = etree.find_event(q)
        #print('event', event, i)

        # perform event and get new events
        old_events, new_events = kmc.move(event, i)

        # update binary search tree
        etree.update_tree(old_events, new_events)

        # if deposition event, check if a new layer is complete
        #if event_i['type'] == 0:
        #    atoms = check_layer()
        #    if len(atoms) > 0:
        #        etree.rebalance_tree(atoms)

        # advance time
        t += -np.log(np.random.random())/Rs

        print('time:', t)

    xyz, box, grain = kmc.get_conf()
    write_cfg('output.xyz', xyz, box, grain)

# end of kmc_pld.py 
