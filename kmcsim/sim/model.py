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
from collections import Counter, defaultdict
#import events
from .events import EventTree

class KMCModel:
    """Class managing kmc moves and event modifications"""

    def __init__(self, latt_type):
        self.latt_type = latt_type
        self.__setup_neighbors()

    def __setup_neighbors(self):
        """Create lists of neighboring (nn & nnn) sites"""

        nbrlist = []

        if self.latt_type == 'fcc':
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
            # NNN
            nbrlist.append(np.array([ 2, 0, 0]))
            nbrlist.append(np.array([-2, 0, 0]))
            nbrlist.append(np.array([ 0, 2, 0]))
            nbrlist.append(np.array([ 0,-2, 0]))
            nbrlist.append(np.array([ 0, 0, 2]))
            nbrlist.append(np.array([ 0, 0,-2]))
        else:
            raise ValueError(f'Chosen {self.latt_type} lattice. Currently only FCC lattice is supported.')

        self.nbrlist = nbrlist


    def make_lattice(self, xyz, box):
        """
        Set up site lables on FCC lattice
        """
    
        # create simulation box (fill with -1 denoting sites on on the FCC lattice
        latt = -1*np.ones((tuple(box)), dtype=int)

        # cycle through lattice sites and mark FCC sites with 0 (empty site)
        for r in product(range(box[0]), range(box[1]), range(box[2])):
            if sum(r) % 2 == 0:
                latt[tuple(r)] = 0
    
        # fill lattice sites with atom ids
        for i, r in enumerate(xyz, start=1):
            assert sum(r) % 2 == 0, f'{r} Atom not on FCC lattice!'
            latt[tuple(r)] = i
    
        self.latt = latt
        self.box = box
        self.xyz = xyz
        self.nat = len(self.xyz)

        # Set grain number of each substrate atom to 0
        self.grain = [0 for _ in range(self.nat)]


    def find_neighbors(self, ri):

        # search nearest neighbors
        # to determine stable sites (needs at least 3)
        neighbors = []
        grain_numbers = []
        for dr in self.nbrlist[0:12]:
            rj = tuple((ri + dr) % self.box)
            iatom = self.latt[rj]
            neighbors.append((iatom, rj))
            if iatom > 0:
                grain_numbers.append(self.grain[iatom-1])

        return neighbors, grain_numbers

    def get_grain(self, grain_numbers):

        grain_counts = Counter(grain_numbers)
        if 0 in grain_counts:
            del grain_counts[0]
        
        if any(grain_counts):
            #g_number = sorted(grain_counts, key=grain_counts.__getitem__)[-1]
            g_number = grain_counts.most_common(1)[0][0]
        else:
            g_number = max(self.grain) + 1

        return g_number


    def init_events(self, rates):
 
        rates = np.array(rates)

        event_list = [[] for _ in range(rates.shape[0])]

        site_dict = defaultdict(list)

        # Deposition event - find vacant sites available for deposition
        for ix, iy in product(range(self.box[0]), range(self.box[1])):

            # find z position
            for iz in range(self.box[2]):
                if self.latt[ix, iy, iz] == 0:
                    ri = np.array([ix, iy, iz], dtype=int)
 
                    # count number of neighbors in the target position
                    neighbors, grain_numbers = self.find_neighbors(ri)

                    # if 3 or more nearest neighbors create a deposition event
                    if len(grain_numbers) > 2:
                        event = {'type':0}
                        event['atom'] = None
                        event['initial'] = None
                        event['final'] = (ix, iy, iz)
                        event_list[0].append(event)

                        # add event information to the site
                        site_dict[(ix, iy, iz)] = [(0, len(event_list[0])-1)]
                    break

        # diffusion events
        for i, ri in enumerate(self.xyz, start=1):
            if ri[2] < 2: continue

            # cycle over neighbor sites
            for dr in self.nbrlist:
                if dr[2] > 0: continue

                rj = (ri + dr) % self.box

                if self.latt[tuple(rj)] == 0:

                    # count number of neighbors in the target position
                    neighbors, grain_numbers = self.find_neighbors(rj)

                    # if 3 or more nearest neighbors create a diffusion event
                    if len(grain_numbers) > 2:
                        event = {'type':1}
                        event['atom'] = i
                        event['initial'] = ri
                        event['final'] = rj
                        event_list[1].append(event)

                        # add event information to the site
                        site_dict[i] = [(1, len(event_list[1])-1)]

        self.event_list =  event_list
        self.site_dict = site_dict

        # Get dictionary of event type counts
        n_events = np.array([len(e) for e in self.event_list])
        print('Number of events:', n_events)

        # Initiate event data structures
        self.etree = EventTree(rates)
        self.etree.update_events(n_events)

        print('events',len(event_list), np.random.choice(np.random.choice(event_list, 2), 5))
        print(self.etree.event_tree[-1])
        print('-------')
        print([[(e['type'], e['final']) for e in et] for et in event_list])
        print('-------')
        print('site dict', site_dict)


    def move(self, event_type, event_number):
        """
        Perform kmc move given by an event (deposition or diffusion)
        and update event list
        """

        # double check if there are some free spaces (just in case - should
        # follow from zero remaining events)
        if len(self.xyz) == self.box[0]*self.box[1]*self.box[2]/2:
            raise ValueError(f'Lattice is full of atoms, no more events possible.')

        event = self.event_list[event_type][event_number]
        old_events = []
        new_events = []
        n_events = self.etree.n_events

        # deposition event
        if event_type == 0:
            ri = event['final']

            # search neighbors and grain numbers
            neighbors, grain_numbers = self.find_neighbors(ri)

            # create a new atom
            self.xyz.append(np.array(ri)) 
            iatom = len(self.xyz)

            # put it on a lattice
            self.latt[tuple(self.xyz[-1])] = len(self.xyz)  # atom number to lattice

            # assign a new grain number to the atom
            self.grain.append(self.get_grain(grain_numbers))


            # Identify old events for removal

            # remove the current deposition event
            old_events.append((event_type, event_number))
            del self.site_dict[ri]
            n_events[0] -= 1

            # remove diffusion events of neighbors
            for i, rj in neighbors:

                # cycle over events associated with the neighbor
                estart = len(old_events)
                for et, en in self.site_dict[i]:
                    if self.event_list[et][en]['final'] == ri:
                        old_events.append((et, en))

                # remove old events from site list
                for e in range(estart, len(old_events)): 
                    self.site_dict[i].remove(old_events[e])
                    n_events[1] -= 1

                if len(self.site_dict[i]) == 0:
                    del self.site_dict[i]

            # Scan potential targets for diffusion and deposition
            for rj in vacancies:

                # create diffusion events for the new atom
                new_events.append((1, iatom, ri, rj))
                n_events[1] += 1

                # count number of neighbors in the target position
                neighbors, grain_numbers = self.find_neighbors(rj)

                # create deposition events around the new atom
                if len(neighbors) > 2:
                    new_events.append((0, rj, None, rj))
                    n_events[0] += 1


        elif event_type == 1: # diffusion
            ri = event['initial']
            rj = event['final']
            iatom = event['atom']

            # search neighbors of the initial state
            neighbors_old, _ = self.find_neighbors(ri)

            # move atom to the new position
            self.latt[tuple(ri)] = 0
            self.latt[tuple(rj)] = iatom+1
            self.xyz[iatom] = rj

            # search neighbors and grain numbers for final state 
            neighbors_new, grain_numbers = self.find_neighbors(rj)

            # all neighbors (old and new)
            neighbors = list(set(neighbors_old + neighbors_new))

            # assign a new grain number to the atom
            self.grain[iatom] = get_grain(grain_numbers)

            # remove the current diffusion event
            old_events.append((event_type, event_number))
            self.site_dict[iatom].remove(old_events[-1])
            if len(self.site_dict[iatom]) == 0:
                del self.site_dict[i]
            n_events[1] -= 1

            # remove site event lists
            for site in neighbors:
                if len(self.site_dict[site[1]]) == 0:
                    del self.site_dict[site[1]]

            # Remove the present event
            old_events.append(event)

            # Change events around the moved atom

            # Add a vacancy-target for old neighbors in the place of ri
            for neighbor in neighbors_old:
                pass

            # Add a vacancy-target for new neighbors in the place of rj
            for neighbor in neighbors_new:
                pass


        # update atom count
        self.nat = len(self.xyz)

        # calculate new n_events

        return n_events


    def get_conf(self):
        return self.xyz, self.box, self.grain


    def advance_time(self):
        """
        Time step of the last event and reset total rates

        Returns
        -------
        dt : float
             Time of the latest event
        """

        dt = -np.log(np.random.random())/self.etree.Rs

        return dt


    def step(self):
        """
        Perform a KMC step.
        """

        # return a random event
        event_type, event_number = self.etree.find_event()

        # perform a step prescribed by the event and return lists of affected events
        n_events = self.move(event_type, event_number)
 
        # update binary search tree
        self.etree.update_events(n_events)
 
