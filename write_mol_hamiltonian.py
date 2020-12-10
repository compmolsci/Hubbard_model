
## Code to write Hamiltonian terms for any given molecule. 
## Author: Valay Agarawal, Dipanjali Halder, Rahul Maitra
## Created 18 August 2020

#Import required stuff, see respective files for definition of each function/variable
import numpy as np
import bk_basic_functions
from multiplying_operators import multiply_two_operators, multiply_four_operators,clean_up
from mol_hamiltonian import onebody_mh, twobody_mh
import inp
import cmath
j=complex(0,1)
from math import floor
from string_interconversion import gen_gates_from_string
transformation = inp.transformation
spin_sites=inp.spin_sites 
V=inp.twoelecint_mo
F=inp.oneelecint_mo
one_body_mh=onebody_mh()
two_body_mh=twobody_mh()
"""
The program is to integrate the PYSCF Integrals and gate structure based on JW/BK encoding, and finally write the full hamiltonian. Output is given in the form of matrix with first column giving the indices and the second column denoting the gate structure corresponding to the coefficients

hamiltonian_two_body_term(term)
Step1 : Term is lifted from a row in the full matrix of Hamiltonian which contains all the relevant combination of the operators.
Step2 : All the operators are multiplied in the gate format
Step3 : Spin sites are converted to spaital sites, and their corresponding integral is taken from the two-electron-integral in MO basis
Step4 : The integral is multiplied by the corresponding coefficients, and the final matrix for one of the two-body ter of the Mol Hamiltonian is returned

hamiltonian_one_body_term(term) 
Step1 : Term is lifted from a row in the full matrix of Hamiltonian which contains all the relevant combination of the operators.
Step2 : All the operators are multiplied in the gate format 
Step3 : Spin sites are converted to spaital sites, and their corresponding integral is taken from the one-electron-integral in MO basis   
Step4 : The integral is multiplied by the corresponding coefficients, and the final matrix for one of the one-body term of the Mol Hamiltonian is returned   

First loop:
The values from the matrix of Hamiltonian terms in the form of operators for one body part are taken iteratively, and stacked. Then, the non-zero and common values are clubbed together, and values less than the threshold are removed.

Second loop
The values from the matrix of Hamiltonian terms in the form of operators for two body part are taken iteratively, and stacked. Then, the non-zero and common values are clubbed together, and values less than the threshold are removed.
"""



#This function will create the two body term
def hamiltonian_twobody_term():
    
    #Intialize the matrix of Hamiltonian 
    tb=np.zeros((0,2))
    #print(two_body_mh)

    #iterate over all the terms of the two body term, see the strucutre of the Hamiiltonian in the file "mol_hamiltonian.py". First term is the ceofficient, and the next four terms are the fermionic operators
    for term in two_body_mh:
        #mUltiply four operators after converting them to gate format
        gates=clean_up(multiply_four_operators(term[1],term[2],term[3],term[4]))
        #Remove zero coefficient terms
        gates=gates[gates[:,0]!=0]
        #print(gates)
        #If there remains a term, 
        if len(gates)!=0:
            #Map Spin orbital to spatial orbital
            p=floor(term[1][0]/2)
            q=floor(term[2][0]/2)
            s=floor(term[3][0]/2)
            r=floor(term[4][0]/2)

            #Get corresponding two body term from the spin-orbital adapted term from pyscf code in inp.py
            coeff=V[p,q,r,s]*inp.total_time

            #If the coefficient is non zero
            if coeff!=0:            
                #Multiply the first columns (that contain the coefficient of the corresponding string of gates
                gates[:,0]=gates[:,0]*coeff*term[0]
            #If the pyscf coeff is zero, terms would be zero
            else:
                gates= [0] 
        else:
            gates= [0]
        #Add them to the hamiltonian
        if len(gates)>1:
            #print(len(gates))
            tb=np.vstack((tb,gates))
    #Add common terms, and return the matrix for two body terms in Molecular Hamiltonian
    #print(tb)
    #print(tb)
    return clean_up((tb))


#Similar procedure to that of one body term
def hamiltonian_onebody_term():
    ob=np.zeros((0,2))
    for term in one_body_mh:
        gates=clean_up(multiply_two_operators(term[0],term[1]))
        gates=gates[gates[:,0]!=0]
        if len(gates)!=0:
            coeff=F[floor(term[0][0]/2),floor(term[1][0]/2)]*inp.total_time
            if coeff!=0:
                gates[:,0]=gates[:,0]*coeff
            else:
                gates= [0]
        else:
            gates=[0]
        if len(gates)>1:
            ob=np.vstack((ob,gates))
    #print(ob)
    return clean_up(ob)

#print(hamiltonian_onebody_term())

#Add both the term of the hamiltonian, and add the common terms,
final_hamilt=clean_up(np.vstack((hamiltonian_onebody_term(),hamiltonian_twobody_term())))
#print(final_hamilt)
#Remove any computational noise
final=final_hamilt[abs(final_hamilt[:,0])>10**(-8)]
#print(final)
def to_dict(final):#,summed_op):
    #print(operator)
    #final=final[0]
    a={}
    key="paulis"
    a.setdefault(key,[])
    for row in final:
        paulis=''
        for element in row[1]:
            paulis+=str(element)
        dictb =  {"coeff": {"imag": row[0].imag, "real": row[0].real }, "label": paulis }
        a[key].append(dictb)
    return a       #paulis.append(element)
        #print(paulis)
#a[key].append(to_dict(final[0]))
#a[key].append(to_dict(final[1]))
summed_op_home=to_dict(final)
#print(summed_op_home['paulis'])
#print(a)
#print(**to_dict(final[0]),**to_dict(final[1]))
#Split the Hamiltnian into two parts that commute within themselves. 
T1=np.array([row for row in final if 'X' in row[1] or 'Y' in row[1]])
T2=np.array([row for row in final if row not in T1])
