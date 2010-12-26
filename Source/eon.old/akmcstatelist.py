##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##
##-----------------------------------------------------------------------------------
""" The statelist module. """

import logging
logger = logging.getLogger('statelist')
import math
import os
import shutil
import sys

from ConfigParser import SafeConfigParser 

import atoms
import config
import akmcstate
import statelist


class AKMCStateList(statelist.StateList):
    """ The StateList class.  Serves as an interface to State objects and StateList metadata. """
    def __init__(self, state_path, kT, thermal_window, max_thermal_window, 
                 epsilon_e, epsilon_r, use_identical, initial_state = None, 
                 list_search_results = False, filter_hole = False):
        statelist.StateList.__init__(self, state_path, epsilon_e, epsilon_r,
                                     use_identical, akmcstate.AKMCState, list_search_results, 
                                     initial_state)
        # aKMC data.
        self.kT = kT
        self.thermal_window = thermal_window
        self.max_thermal_window = max_thermal_window
        self.filter_hole = filter_hole

    def register_process(self, reactant_number, product_number, process_id):
        # Get the reactant and product state objects.
        reactant = self.get_state(reactant_number)
        product = self.get_state(product_number)
        reactant.load_process_table()
        product.load_process_table()
        reverse_procs = product.get_process_table()
        # Make the reactant process point to the product state number.
        reactant.procs[process_id]["product"] = product_number
        reactant.save_process_table()
        if product.get_num_procs() != 0:
            # Check to see if the reverse process is already identified (not -1). If so, return.
            for id in reverse_procs.keys():
                proc = reverse_procs[id]
                if proc["product"] == reactant_number:
                    return
            # The reverse process might already exist, but is unidentified (-1). See if this is the case and fix it.
            candidates = []
            for id in reverse_procs.keys():
                if reverse_procs[id]["product"] == -1:
                    candidates.append(id)
            esaddle = proc["saddle_energy"]
            energetically_close = []
            for id in candidates:
                proc = reverse_procs[id]
                if abs(proc["saddle_energy"] - esaddle) < self.epsilon_e:
                    energetically_close.append(id)
            if len(energetically_close) > 0:
                saddle_config = reactant.get_reactant()
                for id in energetically_close:
                    temp_config = product.get_process_saddle(id)
                    dist = max(atoms.per_atom_norm(saddle_config.r - temp_config.r, saddle_config.box))
                    if dist < self.epsilon_r:
                        reverse_procs[id][product] = product_number
                        product.save_process_table()
                        return
        else:
            # There are no processes, so this must be a new state. Add the reverse process.
            product.set_energy(reactant.procs[process_id]["product_energy"])
            # Update the reactant state to point at the new state id.
            reactant.procs[process_id]['product'] = product_number
            reactant.save_process_table()
            # Create the reverse process in the new state.
            barrier = reactant.procs[process_id]['saddle_energy'] - reactant.procs[process_id]['product_energy']
            shutil.copy(reactant.proc_saddle_path(process_id), product.proc_saddle_path(0))
            shutil.copy(reactant.proc_reactant_path(process_id), product.proc_product_path(0))
            shutil.copy(reactant.proc_product_path(process_id), product.proc_reactant_path(0))
            shutil.copy(reactant.proc_results_path(process_id), product.proc_results_path(0))
            shutil.copy(reactant.proc_mode_path(process_id), product.proc_mode_path(0))
            product.append_process_table(id = 0, 
                                         saddle_energy = reactant.procs[process_id]['saddle_energy'], 
                                         prefactor = reactant.procs[process_id]['product_prefactor'], 
                                         product = reactant_number, 
                                         product_energy = reactant.get_energy(),
                                         product_prefactor = reactant.procs[process_id]['prefactor'],
                                         barrier = barrier, 
                                         rate = reactant.procs[process_id]['product_prefactor'] * math.exp(-barrier / self.kT), 
                                         repeats = 0)
            product.save_process_table()

    def connect_states(self, states):
        for i in states:
            proc_tab = i.get_process_table()
            for j in proc_tab:
                if proc_tab[j]['product'] != -1:
                    continue
                enew = proc_tab[j]['product_energy']
                energetically_close = []
                for state in states:
                    if abs(state.get_energy() - enew) < self.epsilon_e:
                        energetically_close.append(state)

                # Perform distance checks on the energetically close configurations.
                if len(energetically_close) > 0:
                    pnew = i.get_process_product(j)
                    for state in energetically_close:
                        p = state.get_reactant()
                        if self.use_identical:
                            if atoms.identical(p, pnew, self.epsilon_r):
                                # Update the reactant state to point at the new state id.
                                self.register_process(i.number, state.number, j)                            
                        else:
                            dist = max(atoms.per_atom_norm(p.r - pnew.r, p.box))
                            if dist < self.epsilon_r:
                                self.register_process(i.number, state.number, j)                            
