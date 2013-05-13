##-----------------------------------------------------------------------------------
## eOn is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## A copy of the GNU General Public License is available at
## http://www.gnu.org/licenses/
##-----------------------------------------------------------------------------------


import os
import numpy
import logging
logger = logging.getLogger('superbasin')


class Superbasin:
    """Class to manage superbasins: Calculate the mean residence time, exit probabilities, and perform Monte Carlo transitions out of the basin, """\
    """based on Novotny's Absorbing Markov Chain algorithm."""


    def __init__(self, path, id, state_list = None, get_state = None):
        if state_list is None and get_state is None:
            raise ValueError('Superbasin must either have a list of states or a reference to get_state of a StateList')
        self.id = int(id)
        self.path = os.path.join(path, str(self.id))
        if not os.path.isfile(self.path):
            self.states = state_list
            self.state_dict = {}
            for state in state_list:
                self.state_dict[state.number] = state
            self.state_numbers = [state.number for state in state_list]
            self.write_data()
        else:
            self.read_data(get_state)


    def pick_exit_state(self, entry_state):
        """Choose an exit state (state of the basin from which we will be leaving) using absorbing Markov chain theory."""
        for i in range(len(self.state_numbers)):
            if entry_state.number == self.state_numbers[i]:
                entry_state_index = i
                break
        else:
            raise ValueError('Passed entry state is not in this superbasin')

        probability_vector = self.probability_matrix.transpose()[entry_state_index]
        if abs(1.0-numpy.sum(probability_vector)) > 1e-3:
            logger.warning("Probability vector is not 1.0")
            logger.warning('Probability vector ' + str(probability_vector) + " " + str(numpy.sum(probability_vector)))
        probability_vector /= numpy.sum(probability_vector)

        u = numpy.random.random_sample()
        p = 0.0
        for i in range(len(self.states)):
            p += probability_vector[i]
            if p>u:
                exit_state_index = i 
                break
        else:
            logger.warning("Warning: Failed to select exit state; p = " + str(p))
        time = self.mean_residence_times[entry_state_index]
        return time, exit_state_index

    def test(self, entry_state):
        n=len(self.mean_residence_times)
        p=[0]*n
        m=10
        for i in range(m):
            time, exit_state_index=self.pick_exit_state(entry_state)
            p[exit_state_index]+=1
        for i in range(n):
            p[i]=float(p[i])/m
        print 'pick up p:', p
        return time, exit_state_index

    def step(self, entry_state, get_product_state):
    
        # c_i (forming vector c) is the inverse of the sum of the rates for each transient state i
        # Q is the transient matrix of the canonical markov matrix. 
        # R is the recurrent matrix of the canonical markov matrix.
    
        # Build a mapping between transient states and row/column indices. Used in Q and R.
        st2i = {}
        i2st = {}
        index = 0
        for number in self.state_numbers:
            st2i[number] = index
            i2st[index] = number
            index += 1
            
        # Build a mapping between recurrent state identifiers [a (state number, process id) tuple] and 
        # column indices. Used only in R.
        st2col = {}
        col2st = {}
        index = 0
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in procs.iteritems():
                if proc['product'] not in self.state_numbers:
                    st2col[(number, id)] = index
                    col2st[index] = (number, id)
                    index += 1
        
        # Build c.
        c = numpy.zeros(len(self.state_numbers))
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in procs.iteritems():
                c[st2i[number]] += proc['rate']
            c[st2i[number]] = 1.0 / c[st2i[number]]

        # Build Q and R.
        Q = numpy.zeros((len(self.state_numbers), len(self.state_numbers)))
        R = numpy.zeros((len(self.state_numbers), len(col2st)))
        for number in self.state_numbers:
            procs = self.state_dict[number].get_process_table()
            for id, proc in procs.iteritems():
                if proc['product'] in self.state_numbers:
                    Q[st2i[number], st2i[proc['product']]] += proc['rate'] * c[st2i[number]]
                else:
                    R[st2i[number], st2col[(number, id)]] += proc['rate'] * c[st2i[number]]
        
        
        # import pdb; pdb.set_trace()
        
        t, B = self.mcamc(Q, R, c)
        
        b = B[st2i[entry_state.number],:]
        p = 0.0
        u = numpy.random.sample()
        print 'b', b, sum(b)
        print 'u', u
        for i in range(len(b)):
            print 'p', p
            p += b[i]
            if p >= u:
                exit_state_number, exit_proc_id = col2st[i]
                break
        else:
            logger.warning("Warning: Failed to select rate; p = " + str(p))
        
                    
    def mcamc(self, Q, R, c):
        A = numpy.identity(Q.shape[0]) - Q
        t = numpy.linalg.solve(A, c)
        B = numpy.linalg.solve(A, R)
        return t, B


    def contains_state(self, state):
        return state in self.state_dict.values()


    def write_data(self):
        logger.debug('Saving data to %s' % self.path)
        f = open(self.path, 'w')
        for number in self.state_numbers:
            f.write("%d " % number)
        f.close()


    def read_data(self, get_state):
        logger.debug('Reading data from %s' % self.path)
        self.state_numbers = [number for number in open(self.path, 'r').readline().strip().split()]
        self.states = [get_state(number) for number in self.state_numbers]


    def delete(self, storage=None):
        if storage is None:
            logger.debug('Deleting %s' % self.path)
            os.remove(self.path)
        else:
            logger.debug('Storing %s' % self.path)
            path_storage = storage+str(self.id)
            os.rename(self.path, path_storage)
        self.states = None

