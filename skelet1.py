# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 17:37:52 2022

@author: Lin

"""

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
    
    
    #max for DW_2000Q_6
    a_vars = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7']
    b_vars = ['b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7']
    
    """
    a_vars = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8', 'a9', 'a10', 'a11', 'a12', 'a13', 'a14', 'a15', 'a16', 'a17']
    b_vars = ['b0', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'b11', 'b12', 'b13', 'b14', 'b15', 'b16', 'b17']
    """
    
    for lbl in reversed(a_vars):
        a = (a << 1) | sample[lbl]
    for lbl in reversed(b_vars):
        b = (b << 1) | sample[lbl] 
        
    return a,b


P = 21
print("Number to factor is ", P) 


#max for DW_2000Q_6
bP = "{:016b}".format(P)    # "{:06b}" formats for 6-bit binary


"""
#Advance
bP = "{:036b}".format(P)
"""
print("P in binary form ", bP)

print("Dictionary for number P in binary form ")


#max for DW_2000Q_6
p_vars = ['p0', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13', 'p14', 'p15']


"""
#Advance
p_vars = ['p0', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10',
          'p11', 'p12', 'p13', 'p14', 'p15', 'p16', 'p17', 'p18', 'p19',
          'p20', 'p21', 'p22', 'p23', 'p24', 'p25', 'p26', 'p27', 'p28',
          'p29', 'p30', 'p31', 'p32', 'p33', 'p34', 'p35']
# Convert P from decimal to binary
"""


#max for DW_2000Q_6
fixed_variables = dict(zip(reversed(p_vars), "{:016b}".format(P)))

"""
#Advance
fixed_variables = dict(zip(reversed(p_vars), "{:036b}".format(P)))
"""

fixed_variables = {var: int(x) for(var, x) in fixed_variables.items()}

print(fixed_variables)

# Fix product variables

import dwavebinarycsp as dbc


#max for DW_2000Q_6
csp = dbc.factories.multiplication_circuit(8)


"""
#Advance
csp = dbc.factories.multiplication_circuit(18)
"""

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

sampleset = embedding_sampler.sample(bqm, num_reads=1000, label='Factoring 08-41')

#print("Best solution found: \n",sampleset.first.sample)

# To see helper functions, select Jupyter File Explorer View from the Online Learning page
 
a, b = to_base_ten(sampleset.first.sample)

print("Given integer P={}, found factors a={} and b={}".format(P, a, b))

results = response_to_dict(sampleset)
print(results)

# Inspect
dwave.inspector.show(sampleset)