'''
binary_of_all_states(all_states):
bitpacks all the states and returns a array corresponding to all the bitpacked value

check_if_multiple(matrix1, matrix2):
checks if matrix 1 is scalar multiple of matrix2

jw_transfrom(op_def):
takes the operator, gives out one matrix using Z, I and Q/Q_dag

jw_transform(op_def)
takes the operator, gives out two matrices using Z, I, X and Y

apply gates(wavefn, *gates)
applies all the gates passed onto the given wavefn in qubit form, returns wave function in qubit form

write_string_of_gates(op_def):
gives out a string of gates to be multiplied 

gen_string_from_gates(gate_matrix):
converts gate matrix to string


'''





import numpy as np
import cmath
import itertools
import inp
import collections
from bk_basic_functions import stringchar
from basic_functions import classical_to_qubit, qubit_to_classical
j=complex(0,1)
from inp import transformation
#Total number of sites including spins
spatial_sites=inp.spatial_sites
electrons=inp.electrons
spin_sites=2*spatial_sites

#Defining Pauli Gates
X=np.array([[0,1],[1,0]],dtype=complex)
Y=np.array([[0,-j],[j,0]],dtype=complex)
Z=np.array([[1,0],[0,-1]],dtype=complex)
Had=np.array([[1,1],[1,-1]],dtype=complex)/np.sqrt(2)
T=np.array([[1,0],[0,cmath.exp(j*cmath.pi/4)]])
Q=np.add(X,j*Y)/2
Q_dag=np.add(X,-j*Y)/2
I=np.array([[1,0],[0,1]],dtype=complex)

#Arrays for writing the entire JW Transform
A1=np.array(('Z','X','I'),dtype=object)
A2=np.array(('Z','-jY','I'),dtype=object)
A=np.array(('Z','Q','I'),dtype=object)
A_dag1=np.array(('Z','X','I'),dtype=object)
A_dag2=np.array(('Z','jY','I'),dtype=object)
A_dag=np.array(('Z','Qd','I'),dtype=object)

def R(X,theta):
        return np.exp(-j*theta*X/2)






#Change State to single value (Bit packing)
def binary_value(coeff,state): 
    val=1
    for j in range(len(state)):
        val=val+int(state[j])*pow(2,j)
    return val
def binary_of_all_states(all_states):
    return np.array([binary_value(state[0],state[1]) for state in all_states])

#See if given matrix is a scalar multiple of another, 
def check_if_multiple(matrix1, matrix2):
    a=np.matmul(matrix1, np.linalg.inv(matrix2))
    if a[0][0]==a[1][1] and a[0][0]!=0 and a[0][1]==0:
        return a[0][0]
    else:
        return 0


#Code to apply gates on a wavefunction
def apply_gates(coeff,wavefn, gates):
    wavefn=np.matmul(gates, wavefn)
    return coeff,wavefn


#Code to apply an operator on the wavefunction
def apply_operator(coeff,state, operator):                    
    if coeff==0:                                              
        return 0,0
    else:
        
        #Change the state to qubit form
        coeff_0, inp_state=classical_to_qubit(coeff, state)
        
        #change the operator to gate form
        if transformation == 'Jordan-Wigner':
            gate_string1, gate_string2=write_string_of_gates(operator)
        elif transformation == 'Bravyi-Kitaev':
            gate_string1, gate_string2=stringchar(operator)
        
        #Generate matrix from the gate information
        coeff1,gate_matrix1=gen_gates_from_string(gate_string1)
        coeff2,gate_matrix2=gen_gates_from_string(gate_string2)
        
        #Apply gates
        c1,appl1=apply_gates(coeff,inp_state,gate_matrix1)
        c2,appl2=apply_gates(coeff,inp_state,gate_matrix2)
        
        #Change back to Fermionic form
        a1,s1=qubit_to_classical(c1,appl1)
        a2,s2=qubit_to_classical(c2,appl2)
        return(a1*coeff1, s1, a2*coeff2, s2) 
def apply_operators(coeff,state, *operators): 
    for operator in operators:                                  ####
        coeff,state=apply_operator(coeff, state, operator)      ####
    return coeff, state                                         ####
def apply_string(coeff,state, string):
    if coeff==0:
        return 0,0
    else:
        coeff_0, inp_state=classical_to_qubit(coeff,state)
        coeff1,gate_matrix1=gen_gates_from_string(string)
        c1, appl1=apply_gates(coeff_0, inp_state, gate_matrix1)
        a1,s1=qubit_to_classical(c1, appl1)
        return (a1*coeff1, s1)
##### DEPRECATED


#Return string of gates for a particular creation and annihilation operator
def write_string_of_gates(op_def):
    #print(op_def)
    if inp.arrange=='abab':
        if op_def[0]<spatial_sites:
            k=2*op_def[0]+1
        else:
            k=2*(op_def[0]-spatial_sites)
    else:
    #The operator has two parts of information, the first part defines the target site number, and the second type gives the information about the presence of dagger
        k=op_def[0]
    dag=op_def[1]
    #print(op_def,k,dag)
    #initialize two strings
    string1=[]
    string2=[]


    for l in range(spin_sites):
        #Apply a series of Z gates before the target site
        if l<k:
            string1.append('Z')
            string2.append('Z')

        #For target site use Q, Qdagger
        elif l==k:
            if dag==True:
                string1.append('X')
                string2.append('Y')
                a=-0.5*j
            else:
                string1.append('X')
                string2.append('Y')
                a=0.5*j

        #Rest sites can be just I
        else:
            string1.append('I')
            string2.append('I')
    return((0.5,string1),(a,string2))


#Take the string of gates in list form and create a hamiltonian matrix
def gen_gates_from_string(string):
    coeff=string[0]
    string=string[1]
    #Read the values in the string, and change them into the matrix  form, by a simple one to one mapping
    matrix=[]
    for item in string:
        if item=='Z':
            matrix.append(Z)
        elif item=='X':
            matrix.append(X)
        elif item=='I':
            matrix.append(I) 
        elif item=='Y':
            matrix.append(Y)
    return (coeff,np.array(matrix))

## Change a gate matrix into string of gates in terms of X,y,Z,I
def gen_string_from_gates(gate_matrix):
    string=[]
    sign=1
    coeff=gate_matrix[0]

    #Scan over matrices, and if the matrix is a scalar multiple of the standard Pauli Gates, append it in a string, and the coefficient is multiplied in the coefficient/ 
    for matrix1 in gate_matrix[1]:
        if check_if_multiple(matrix1,X): 
            string.append('X')
            sign=sign*check_if_multiple(matrix1, X)
        elif check_if_multiple(matrix1,Z): 
            string.append('Z')
            sign=sign*check_if_multiple(matrix1,Z)
        elif check_if_multiple(matrix1,I):
            string.append('I')
            sign=sign*check_if_multiple(matrix1,I)
        elif check_if_multiple(matrix1,Y):
            string.append('Y')
            sign=sign*check_if_multiple(matrix1, Y)
        else: 
            print('i dont know')
            print(matrix1)
    return sign*coeff, string
