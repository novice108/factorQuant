#    Copyright 2018 D-Wave Systems Inc.

#    Licensed under the Apache License, Version 2.0 (the "License")
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at

#        http: // www.apache.org/licenses/LICENSE-2.0

#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

import sys
import time
import logging
import functools
from collections import OrderedDict

import dwavebinarycsp as dbc
from dwave.system import DWaveSampler, EmbeddingComposite
#import dwave.inspector

log = logging.getLogger(__name__)

def sanitised_input(description, variable, range_):
    start = range_[0]
    stop = range_[-1]

    while True:
        ui = input("Input {:15}({:2} <= {:1} <= {:2}): ".format(description, start, variable, stop))

        try:
            ui = int(ui)
        except ValueError:
            print("Input type must be int")
            continue

        if ui not in range_:
            print("Input must be between {} and {}".format(start, stop))
            continue

        return ui

def validate_input(ui, range_):
    start = range_[0]
    stop = range_[-1]

    if not isinstance(ui, int):
        raise ValueError("Input type must be int")

    if ui not in range_:
        raise ValueError("Input must be between {} and {}".format(start, stop))

def factor(P):

    # Construct circuit
    # =================
    construction_start_time = time.time()

    validate_input(P, range(2 ** 6))

    # Constraint satisfaction problem
    csp = dbc.factories.multiplication_circuit(3)

    # Binary quadratic model
    bqm = dbc.stitch(csp, min_classical_gap=.1)

    # multiplication_circuit() creates these variables
    p_vars = ['p0', 'p1', 'p2', 'p3', 'p4', 'p5']

    # Convert P from decimal to binary
    fixed_variables = dict(zip(reversed(p_vars), "{:06b}".format(P)))
    fixed_variables = {var: int(x) for(var, x) in fixed_variables.items()}

    # Fix product qubits
    for var, value in fixed_variables.items():
        bqm.fix_variable(var, value)

    log.debug('bqm construction time: %s', time.time() - construction_start_time)

    # Run problem
    # ===========

    sample_time = time.time()
    
    """
    num_reads = 100
    from dwave.system import LeapHybridSampler
    sampleset = LeapHybridSampler().sample(bqm, label='tk')
    print(sampleset)
    """
    
    # Set a QPU sampler
    sampler = EmbeddingComposite(DWaveSampler())

    num_reads = 1000 #Number of runs on QC
    chainstrength = 10 # update
    
    #annealing_time=1000 #var1 - annealing time - Block
    
    """
    # Pause
    anneal_time = 10.0
    pause_duration = 500.0      # Must be greater than 0
    pause_start = 0.4        # Must be between 0 and 1
    schedule=[[0.0,0.0],[pause_start*anneal_time,pause_start],[pause_start*anneal_time+pause_duration, pause_start],[anneal_time+pause_duration, 1.0]]
    # End Pause block
    """
    
    """
    # Quench
    anneal_time = 50.0
    quench_slope = 1.0      # Must be greater than 0
    quench_start = 0.2      # Must be between 0 and 1
    schedule=[[0.0,0.0],[quench_start*anneal_time,quench_start],[(1-quench_start+quench_slope*quench_start*anneal_time)/quench_slope, 1.0]]
    # End of Quench Block
    """
    import numpy as np
    get_hamming_distance = lambda x1, x2: np.sum(x1 != x2)

    def get_hamming_distances(sols):
        sols = np.array(sols)
        return np.array([get_hamming_distance(x1, x2) for x1, x2 in zip(sols, sols[1:])])
    
    from dimod import ising_energy
    import random
    
    (h, J, o) = bqm.to_ising()
    print("h=",h)
    print("J=",J)
    print("o=",o)
    lenh = len(h)
    #print("lenh =", lenh)
    ground_state = {s_i:-1 for s_i in range(lenh)} # Global minimum is for all qubits with spin value -1  
    print("ground_state = ", ground_state)
    print("qubo = ", bqm.to_qubo())
    #ground_energy = ising_energy(ground_state, h, J, o)
    #print("ground_energy = ", ground_energy)

    #print("Energy of the global minimum: {}".format(ground_energy))
    #print("Energy of a random sample: {}".format(ising_energy({s_i: 2*random.randint(0, 1)-1 for s_i in range(lenh)}, h, J)))
    
        
    
    def analyze(answer):
            reads=num_reads
            solutions, energies = answer.record.sample, answer.record.energy
            energy_best = round(answer.first.energy, 2)
            #ratio = list(answer.record.energy).count(ground_energy)/float(reads)     
            hamming_distances = get_hamming_distances(solutions)
            energy_mean = round(energies.mean(), 2)
            hamming_mean = round(hamming_distances.mean(), 2)
        #return([solutions, energies, hamming_distances, energy_best, ratio, energy_mean, hamming_mean])
            return([solutions, energies, hamming_distances, energy_best, energy_mean, hamming_mean])

    
    sampleset = sampler.sample(bqm,
                               num_reads=num_reads,
                               #anneal_schedule=schedule,
                               chain_strength = chainstrength,
                               #annealing_time=1000 - Block
                               label='Schedule + Factoring')
    #print("analize ", analyze(sampleset))
    #dwave.inspector.show(sampleset)
    log.debug('embedding and sampling time: %s', time.time() - sample_time)


    # Output results
    # ==============
    output = {
        "Results": [],
     #   {
     #       "a": Number,
     #       "b": Number,
     #       "Valid": Boolean,
     #       "Occurrences": Number,
     #       "Percentage of results": Number
     #    }
    "Timing": {
        "Actual": {
            "QPU processing time": None  # microseconds
        }
    },
    "Number of reads": None
    }

    # multiplication_circuit() creates these variables
    a_vars = ['a0', 'a1', 'a2']
    b_vars = ['b0', 'b1', 'b2']

    results_dict = OrderedDict()
    for sample, num_occurrences in sampleset.data(['sample', 'num_occurrences']):
        # Convert A and B from binary to decimal
        a = b = 0
        for lbl in reversed(a_vars):
            a = (a << 1) | sample[lbl]
        for lbl in reversed(b_vars):
            b = (b << 1) | sample[lbl]
        # Cast from numpy.int to int
        a, b = int(a), int(b)
        # Aggregate results by unique A and B values (ignoring internal circuit variables)
        if (a, b, P) in results_dict:
            results_dict[(a, b, P)]["Occurrences"] += num_occurrences
            results_dict[(a, b, P)]["Percentage of results"] = 100 * \
                results_dict[(a, b, P)]["Occurrences"] / num_reads
        else:
            results_dict[(a, b, P)] = {"a": a,
                                       "b": b,
                                       "Valid": a * b == P,
                                       "Occurrences": num_occurrences,
                                       "Percentage of results": 100 * num_occurrences / num_reads}

    output['Results'] = list(results_dict.values())
    output['Number of reads'] = num_reads

    output['Timing']['Actual']['QPU processing time'] = sampleset.info['timing']['qpu_access_time']

    return output

def display_output(output):
    header1_str = 'Factors    Valid?  Percentage of Occurrences'
    header2_str = ' ' * header1_str.index('P') + 'Numeric & Graphic Representation'
    total_width = 80  # Assumed total console width
    # Width available to draw bars:
    available_width = total_width - header1_str.index('P') - 4

    header_len = max(len(header1_str), len(header2_str))
    print('-'*header_len)
    print(header1_str)
    print(header2_str)
    print('-'*header_len)

    for result in output['Results']:
        percentage = result['Percentage of results']
        print('({:3},{:3})  {:3}     {:3.0f} '.format(result['a'], result['b'], 'Yes' if result['Valid'] else '', percentage), end='')
        nbars = int(percentage/100 * available_width)
        print('*' * nbars)

if __name__ == '__main__':
    # get input from user
    print("Enter a number to be factored:")
    P = sanitised_input("product", "P", range(2 ** 6))

    # send problem to QPU
    print("Running on QPU")
    output = factor(P)

    # output results
    #display_output(output)
    