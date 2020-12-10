import numpy as np
import inp
import itertools
spatial_sites=inp.spatial_sites
spin_sites=2*spatial_sites
electrons=inp.electrons
j=complex(0,1)
occ=np.array([[0],[1]])
unocc=np.array([[1],[0]])
'''

Has the following functions

generate_possible_states():  
Takes no input, generate all states based on number of sites and electrons

classical_to_qubit(string_of_occupancy_numbers):
converts occupancy number string into a matrix of qubits

gen_all_operators():
Takes no input, generates all operators in the form of a operator number and dagger identifier

qubit_to_classical(matrix_of_qubits):
Converts matrix of qubits into a occupancy number string
'''


def generate_possible_states():
    all_states=np.array(list(itertools.product([1,0],repeat=spin_sites)))
    states1=np.array([np.sum(state) for state in all_states])
    all_states=all_states[states1==electrons]
    all_States=[]
    for state in all_states:
        all_States.append(np.array((1,state)))
    return np.array(all_States)

def classical_to_qubit(coeff,string):
    '''
    This functions inputs a string of occupation in fermionic basis, along with the coefficient of the string and changes it into qubit basis. It returns a matrix of the form 1*2*n, where n is the number of qubits(which is equal to the spin sites)
    '''
    #Initialize the matrix 
    occ_matrix=[]
    
    #scan along a string
    for char in string: 
        #If character in string is 0, then the site is empty, and we have to initialize it with |0>
        if char==0:
            occ_matrix.append(unocc)

        #If character in string is 1, then the site if filled, and we have to initialize it with |1>
        else:
            occ_matrix.append(occ)
    return coeff,np.array(occ_matrix)

def gen_all_operators():
    s=spin_sites
    values=[]
    for site in range(s):
        values.append(np.array((site,1)))
        values.append(np.array((site,0)))
    return np.array(values)

def qubit_to_classical(coeff,occ_matrix):
    '''
    This function inputs a matrix of occupation in qubit basis, and changes it back to the fermionc basis. However, this function is now depracated'''

    sign=1*coeff
    #Initialize a string of Fermionic occupation 
    string =[]
    
    #If the coefficient of the matrix is 0, then there is no state, and return 0
    if coeff==0:
        return 0,0

    #Scan along the matrix, for each 1*2 matrix, and if the matrix is of the form [0,1], append 1 in the string, else if it is of the form [1,0], append 1
    else:
        for site in occ_matrix:
            if site[0]==0 and site[1]!=0:
                string.append('1')
                sign=sign*site[1]
            elif site[1]==0 and site[0]!=0:
                string.append(0)
                sign=sign*site[0]
            else: 
                1
                print(site)
        return sign,string
