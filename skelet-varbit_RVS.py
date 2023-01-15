# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 17:37:52 2022

@author: Lin

"""
import sys
import time
import logging
import functools
from collections import OrderedDict

def response_to_dict(sampleset, M):
    results_dict = OrderedDict()
    for sample, energy in sampleset.data(['sample', 'energy']):
        a, b = to_base_ten(sample, M)
        if (a, b) not in results_dict:
            results_dict[(a, b)] = round(energy, 2)
            
    return results_dict

def to_base_ten(sample, M):
    a = b = 0
    
    # we know that multiplication_circuit() has created these variables
    
    # Variables for output factors a and b
    a_vars = dict()
    for i in range(int(M/2)):
        a_vars['a'+str(i)] = i
    #print(a_vars)
    
    b_vars = dict()
    for i in range(int(M/2)):
        b_vars['b'+str(i)] = i
    #print(b_vars)
    
    for lbl in reversed(a_vars):
        a = (a << 1) | sample[lbl]
    for lbl in reversed(b_vars):
        b = (b << 1) | sample[lbl] 
        
    return a,b

log = logging.getLogger(__name__)
sample_time = time.time()
log.debug('embedding and sampling time: %s', time.time() - sample_time)

# Output results
# ==============

output = {
    "Results": [],
    #    {
    #        "a": Number,
    #        "b": Number,
    #        "Valid": Boolean,
    #        "Occurrences": Number,
    #        "Percentage of results": Number,
    #        "Temperature": Number
    #    }
    "Timing": {
        "Actual": {
            "QPU processing time": None  # microseconds
        }
    },
    "Number of reads": None
}

Pa = 3 
Pb = 7
import math

P = Pa*Pb # Input number P for factoring in decimal format
print("Number to factor is ", P, "=", Pa, '*', Pb)
#print("Number to factor in binary", bin(P))


a = 0
M = 0

if Pa >= Pb :
    g = math.log2(Pa)
    #print(g)
    a = int(g)
    if (g%1) > 0:
        a = int(g) + 1
    print(a)
    if (a%2) != 0:
        a = a + 1
    print(a)

if Pb > Pa :
    g = math.log2(Pb)
    print(g)
    a = int(g)
    if (g%1) > 0:
        a = int(g) + 1
    print(a)
    if (a%2) != 0:
        a = a + 1
    print(a)

M = 2*a
if M < 6:
    M = 6

print("M", M)

  
z = '{'+':'+'0'+str(M)+'b'+'}'
#print(z)
#print(z.format(P))
p_vars = {}
for i in range(M):
    p_vars['p'+str(i)] = bin(0)
#print(p_vars)


# input number P in binary format
#bP = "{:06b}".format(P)    # "{:06b}" formats for 6-bit binary
#print("P in binary form ", bP)


# Convert P from decimal to binary and put it in dictionary fixed_variables
fixed_variables = dict(zip(reversed(p_vars), z.format(P)))
#print("fv", fixed_variables)
fixed_variables = {var: int(x) for(var, x) in fixed_variables.items()}
#print("Dictionary for number P in binary form ")
#print(fixed_variables)


# Fix product variables
import dwavebinarycsp as dbc
csp = dbc.factories.multiplication_circuit(int(M/2))
bqm = dbc.stitch(csp, min_classical_gap=.1)

for var, value in fixed_variables.items():
    bqm.fix_variable(var, value)


#print("BQM has {} non-fixed variables: \n\t{}".format(len(bqm.variables), list(bqm.variables)))


from dwave.system import DWaveSampler
import dwave.inspector

# Use a D-Wave system as the sampler
sampler = DWaveSampler() 

print("QPU {} was selected.".format(sampler.solver.name))

from dwave.system import EmbeddingComposite

embedding_sampler = EmbeddingComposite(sampler)


from datetime import datetime

now = datetime.now() # current date and time

date_time = now.strftime("%m/%d/, %H:%M")
print("date and time:",date_time)	

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

num_reads = 100 #Number of runs on QC
chainstrength = 10 # update

max_slope = 1.0/DWaveSampler().properties["annealing_time_range"][0]
    
reverse_schedule = make_reverse_anneal_schedule(s_target=0.45, hold_time=80, ramp_up_slope=max_slope)
time_total = reverse_schedule[3][0]

from dwave.system import DWaveSampler, EmbeddingComposite
sampler = EmbeddingComposite(DWaveSampler())

forward_answer = sampler.sample(bqm,num_reads=num_reads,
                                annealing_time=time_total,
                                chain_strength = chainstrength,
                                label='Forward Reverse Annealing' + str(date_time) + " " + str(P)
                                     + " Factoring" + str(num_reads)+" runs",
                                     answer_mode='raw')
forward_solutions, forward_energies = forward_answer.record.sample, forward_answer.record.energy
i5 = int(5.0/95*len(forward_answer))  # Index i5 is about a 5% indexial move from the sample of lowest energy

initial = dict(zip(forward_answer.variables, forward_answer.record[i5].sample))

reverse_anneal_params = dict(anneal_schedule=reverse_schedule, initial_state=initial, reinitialize_state=False) #True
num_reads = num_reads*10
reverse_answer = sampler.sample(bqm, 
                                  num_reads=num_reads, 
                                  label='Reverse Annealing'+str(date_time) + " " + str(P) +'='+str(Pa)+'*'+str(Pb)+ " Factoring" + str(num_reads)+" runs",
                                  answer_mode='raw', #'histogram' 
                                  **reverse_anneal_params)
reverse_solutions, reverse_energies = reverse_answer.record.sample, reverse_answer.record.energy

#print("reverse_energies ", reverse_energies)

sampleset = reverse_answer
print("sampleset")
print(sampleset)

"""
#without reverse annealing
from dwave.system import EmbeddingComposite
embedding_sampler = EmbeddingComposite(sampler)
sampleset = sampler.sample(bqm, num_reads=num_reads, label = str(date_time) + " " + str(P)+'='+str(Pa)+'*'+str(Pb))
                                     + " Factoring" + str(num_reads)+" runs")
"""

#print("Best solution found: \n",sampleset.first.sample)

# To see helper functions, select Jupyter File Explorer View from the Online Learning page
 
a, b = to_base_ten(sampleset.first.sample, M)

print("Given integer P={}, found factors a={} and b={}".format(P, a, b))

results = response_to_dict(sampleset, M)
#print(results)

results_dict = OrderedDict()

a_vars = dict()
for i in range(int(M/2)):
    a_vars['a'+str(i)] = i
#print(a_vars)

b_vars = dict()
for i in range(int(M/2)):
    b_vars['b'+str(i)] = i
#print(b_vars)


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

def display_output(output, results):
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

#https://stackoverflow.com/questions/15114843/accessing-dictionary-value-by-index-in-python
    m = 0
    for result in output['Results']:    
        percentage = result['Percentage of results']
        print('({:3},{:3})  {:3}     {:3.0f} '.format(result['a'], result['b'], 'Yes' if result['Valid'] else '', percentage), end='')
        nbars = int(percentage/100 * available_width)
        print('*' * nbars), print(list(results.values())[m])
        m = m + 1

display_output(output, results)

# Inspect
#dwave.inspector.show(sampleset.aggregate()) #did not convert from raw to histogram
#dwave.inspector.show(sampleset)
dwave.inspector.show(forward_answer, block='never')



"""
bidict = {}
N = 10
for i in range(N):
    bidict['p'+str(i)] = i
print(bidict)

M = 6
z = '{'+':'+'0'+str(M)+'b'+'}'
print(z)
print(z.format(21))

import math

g = math.log2(9)
print(g)
if (g%1) > 0:
    g = g + 1
print(int(g))
"""