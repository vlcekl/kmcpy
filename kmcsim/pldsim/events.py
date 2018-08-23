#!//anaconda/envs/py36/bin/python
#
# File name:   kmc_pld.py
# Date:        2018/08/03 09:07
# Author:      Lukas Vlcek
#
# Description: 
#

import numpy as np

class EventTree:
    """
    Class maintaining a binary tree for random event type lookup
    and arrays for choosing specific event. 
    """

    def __init__(self, rates):

        # list of rates
        self.rates = rates

        # empty event type tree
        self.t = []

        # empty event id (list of lists)
        self.eid = [[] for _ in range(len(rates))]


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
        return self.t[-1]


    def update_events(self, old_events, new_events):
        """
        Update tree: remove old events and add new events
        """
        pass


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
