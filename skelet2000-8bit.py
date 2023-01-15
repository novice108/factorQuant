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

def response_to_dict(sampleset):
    results_dict = OrderedDict()
    for sample, energy in sampleset.data(['sample', 'energy']):
        a, b = to_base_ten(sample)
        if (a, b) not in results_dict:
            results_dict[(a, b)] = round(energy, 2)
            
    return results_dict

def to_base_ten(sample):
    a = b = 0
    
    # we know that multiplication_circuit() has created these variables
    
    # Variables for output factors a and b
    a_vars = ['a0', 'a1', 'a2', 'a3']
    b_vars = ['b0', 'b1', 'b2', 'b3']

   
    
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

P = 10 # 55 = 11x5 input number P for factoring in decimal format
print("Number to factor is ", P) 


# input number P in binary format
bP = "{:08b}".format(P)    # "{:06b}" formats for 6-bit binary


print("P in binary form ", bP)

print("Dictionary for number P in binary form ")
p_vars = ['p0', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']

# Convert P from decimal to binary and put it in dictionary fixed_variables
fixed_variables = dict(zip(reversed(p_vars), "{:08b}".format(P)))
fixed_variables = {var: int(x) for(var, x) in fixed_variables.items()}
print(fixed_variables)

"""
# Fix product variables
import dwavebinarycsp as dbc
csp = dbc.factories.multiplication_circuit(3)
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

num_reads = 100
sampleset = embedding_sampler.sample(bqm, num_reads=100, label='2x5 Factoring 100 runs')

#print("Best solution found: \n",sampleset.first.sample)

# To see helper functions, select Jupyter File Explorer View from the Online Learning page
 
a, b = to_base_ten(sampleset.first.sample)

print("Given integer P={}, found factors a={} and b={}".format(P, a, b))

results = response_to_dict(sampleset)
print(results)

results_dict = OrderedDict()
a_vars = ['a0', 'a1', 'a2']
b_vars = ['b0', 'b1', 'b2']
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
#dwave.inspector.show(sampleset)
"""


"""
bidict = {}
N = 10
for i in range(N):
    bidict['p'+str(i)] = i
print(bidict)
"""

"""
M = 6
z = '{'+':'+'0'+str(M)+'b'+'}'
print(z)
print(z.format(21))
"""

"""
import math

g = math.log2(9)
print(g)
if (g%1) > 0:
    g = g + 1
print(int(g))
"""