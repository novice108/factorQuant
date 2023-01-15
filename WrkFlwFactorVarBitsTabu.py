# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 05:41:03 2022
@author: Lin
Trial hybrid with various lenght of number to find factors
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

Pa = 311 
Pb = 313
import math

P = Pa*Pb # Input number P for factoring in decimal format
print("Number to factor is "+ str(P)+'='+str(Pa)+'*'+str(Pb))
print("Number to factor in binary", bin(P))


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

# Convert P from decimal to binary and put it in dictionary fixed_variables
fixed_variables = dict(zip(reversed(p_vars), z.format(P)))
#print("fv", fixed_variables)
fixed_variables = {var: int(x) for(var, x) in fixed_variables.items()}
print("Dictionary for number P in binary form ")
print(fixed_variables)


# Fix product variables
import dwavebinarycsp as dbc
csp = dbc.factories.multiplication_circuit(int(M/2))
bqm = dbc.stitch(csp, min_classical_gap=.1)

for var, value in fixed_variables.items():
    bqm.fix_variable(var, value)


#print("BQM has {} non-fixed variables: \n\t{}".format(len(bqm.variables), list(bqm.variables)))
from datetime import datetime
now = datetime.now() # current date and time
date_time = now.strftime("%m/%d/, %H:%M")
print("date and time:",date_time)	
#from dwave.system import LeapHybridSampler
import hybrid
state = hybrid.State.from_problem(bqm)
workflow = (hybrid.Parallel(hybrid.TabuProblemSampler(num_reads=10000, timeout=50),
        hybrid.SimulatedAnnealingProblemSampler(num_reads=10000)) |
        hybrid.ArgMin() )
    
sampleset1 = workflow.run(state).result()
#print(sampleset1)
sampleset = sampleset1['samples']
#print("wsamplset", sampleset)
#sampleset = LeapHybridSampler().sample(bqm, label="LeapHybrid "+str(date_time)+" "+ str(P)+str('=')+str(Pa)+str('*')+str(Pb)+" Factoring")
#print("Lsampleset", sampleset)

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

num_reads = 100
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
#output['Number of reads'] = num_reads

#output['Timing']['Actual']['QPU processing time'] = sampleset.info['timing']['qpu_access_time']

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
