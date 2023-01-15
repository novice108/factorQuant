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

#sampleset1 = 1

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

def make_reverse_anneal_schedule(s_target=0.0, hold_time=10.0, ramp_back_slope=0.2, ramp_up_time=0.0201,
                                 ramp_up_slope=None):
    """Build annealing waveform pattern for reverse anneal feature.

    Waveform starts and ends at s=1.0, descending to a constant value
    s_target in between, following a linear ramp.

      s_target:   s-parameter to descend to (between 0 and 1)
      hold_time:  amount of time (in us) to spend at s_target (must be >= 2.0us)
      ramp_slope: slope of transition region, in units 1/us
    """
    # validate parameters
    if s_target < 0.0 or s_target > 1.0:
        raise ValueError("s_target must be between 0 and 1")
    if hold_time < 0.0:
        raise ValueError("hold_time must be >= 0")
    if ramp_back_slope > 0.2:
        raise ValueError("ramp_back_slope must be <= 0.2")
    if ramp_back_slope <= 0.0:
        raise ValueError("ramp_back_slope must be > 0")

    ramp_time = (1.0 - s_target) / ramp_back_slope

    initial_s = 1.0
    pattern = [[0.0, initial_s]]

    # don't add new points if s_target == 1.0
    if s_target < 1.0:
        pattern.append([round(ramp_time, 4), round(s_target, 4)])
        if hold_time != 0:
            pattern.append([round(ramp_time+hold_time, 4), round(s_target, 4)])

    # add last point
    if ramp_up_slope is not None:
        ramp_up_time = (1.0-s_target)/ramp_up_slope
        pattern.append([round(ramp_time + hold_time + ramp_up_time, 4), round(1.0, 4)])
    else:
        pattern.append([round(ramp_time + hold_time + ramp_up_time, 4), round(1.0, 4)])

    return pattern


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
    #sampler = EmbeddingComposite(DWaveSampler())
    #sampler = EmbeddingComposite()
    from dwave.system import DWaveSampler
    sampler = DWaveSampler()
    
    num_reads = 1000 #Number of runs on QC
    chainstrength = 20 # update
    
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
    
    max_slope = 1.0/DWaveSampler().properties["annealing_time_range"][0]
    
    reverse_schedule = make_reverse_anneal_schedule(s_target=0.45, hold_time=80, ramp_up_slope=max_slope)
    time_total = reverse_schedule[3][0]
    
    from dwave.system import DWaveSampler, EmbeddingComposite
    sampler = EmbeddingComposite(DWaveSampler())
    #print("QPU {} was selected.".format(sampler.solver.name))
    
    forward_answer = sampler.sample(bqm,num_reads=1000,annealing_time=time_total,label='Forward Annealing 1000',answer_mode='raw')
    forward_solutions, forward_energies = forward_answer.record.sample, forward_answer.record.energy
    i5 = int(5.0/95*len(forward_answer))  # Index i5 is about a 5% indexial move from the sample of lowest energy

    initial = dict(zip(forward_answer.variables, forward_answer.record[i5].sample))

    reverse_anneal_params = dict(anneal_schedule=reverse_schedule, initial_state=initial, reinitialize_state=False)
    reverse_answer = sampler.sample(bqm, 
                                      num_reads=1000, 
                                      label='Reverse Annealing 1000',
                                      answer_mode='raw', 
                                      **reverse_anneal_params)
    reverse_solutions, reverse_energies = reverse_answer.record.sample, reverse_answer.record.energy
     
    sampleset = reverse_answer
    #sampleset1 = sampleset
    #print(sampleset)
    
    """
    sampleset = sampler.sample(bqm,
                               num_reads=num_reads,
                               #anneal_schedule=schedule,
                               chain_strength = chainstrength,
                               #annealing_time=1000 - Block
                               label='Schedule + Factoring')
    """
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

    num_reads = 1000
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
    
    """
    # Inspect
    import dwave.inspector
    dwave.inspector.show(sampleset1)
    """

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
    display_output(output)
    
