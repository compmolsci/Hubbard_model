### Code to multiply operators and gather common terms
### Author: Valay Agarawal, Dipanjali Halder
### Date: 18 Aug 2020



import numpy as np
import inp
import cmath
spatial_sites=inp.spatial_sites
from string_interconversion import write_string_of_gates, gen_gates_from_string, gen_string_from_gates,apply_gates
from basic_functions import qubit_to_classical, classical_to_qubit
from bk_basic_functions import stringchar
j=complex(0,1)
'''
multiplying_two_operators(operator1, operator2):
    Takes two operator and returns a matrix in the form of gates containing all the terms(4 terms)
    Step1 : See the appropriate transformation
    Step2 : Convert operator from operator form to string of gates
    Step3 : Multiply coefficients and gates
    Step4 : Convert matrix of gates into string of gates

multiplying_four_operators(operator1, operator2, operator3, operator4):
    Takes two operator and returns a matrix in the form of gates containing all the terms(16 terms)
    Same as multiplying_two_operators

clean_up_multiplication(gates1, gates2):
    Takes two matrices of gates and reduce them to addition form
    
clean_up(gates):
    Takes one matrix of gates and reduce to only the unique gate strings and add common terms
'''
transformation=inp.transformation
def multiply_two_operators(operator1, operator2):
    product=[]
    if transformation=='Jordan-Wigner':
        string1=write_string_of_gates(operator1)
        string2=write_string_of_gates(operator2)
    elif transformation=='Bravyi-Kitaev':
        string1=stringchar(operator1)
        string2=stringchar(operator2)
    for gates1 in string1:
        for gates2 in string2:     
            A1=gen_gates_from_string(gates1)
            A2=gen_gates_from_string(gates2)
            product.append(((A1[0]*A2[0]),np.matmul(A1[1],A2[1])))
    gates=[]
    for gate_matrix in product:
        gates.append(gen_string_from_gates(gate_matrix))   
    return np.array(gates)

def multiply_four_operators(operator1,operator2,operator3, operator4):
    product=[]
    if transformation=='Jordan-Wigner': 
        string1=write_string_of_gates(operator1)
        string2=write_string_of_gates(operator2)
        string3=write_string_of_gates(operator3)
        string4=write_string_of_gates(operator4)
    elif transformation=='Bravyi-Kitaev':
        string1=stringchar(operator1)
        string2=stringchar(operator2)
        string3=stringchar(operator3)
        string4=stringchar(operator4)
    for gates1 in string1:
        for gates2 in string2:
            for gates3 in string3:
                for gates4 in string4:
                    A1=gen_gates_from_string(gates1)
                    A2=gen_gates_from_string(gates2)
                    A3=gen_gates_from_string(gates3)
                    A4=gen_gates_from_string(gates4)
                    product.append(((A1[0]*A2[0]*A3[0]*A4[0]),np.matmul(np.matmul(np.matmul(A1[1],A2[1]),A3[1]),A4[1])))
    gates=[]
    for gate_matrix in product:
        gates.append(gen_string_from_gates(gate_matrix))
    #print(np.array(gates))
    return np.array(gates)

def clean_up_multiplication(gates1, gates2):
    final=[]
    for gate_string1 in gates1:
        for gate_string2 in gates2:
            if np.array_equal(np.array(gate_string1[1]),np.array(gate_string2[1])) and gate_string1[0]-gate_string2[0]!=0:
                final.append(np.array((gate_string1[0]-gate_string2[0],gate_string1[1])))
    return np.array(final)

def clean_up(gates):
    gate_unique=np.unique(gates[:,1])
    t=np.zeros((len(gate_unique)),dtype=complex)
    for i in range(len(gate_unique)):
        for j in range(len(gates)):
            if(np.array_equal(np.array(gates[j,1]),np.array(gate_unique[i]))):
                t[i]+=gates[j,0]
    t=np.transpose(np.vstack((t,gate_unique)))
    return t[t[:,0]!=0]
