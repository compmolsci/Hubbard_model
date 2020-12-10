
### Code to write Hubbard Hamiltonian based on JW/BK encoding in 1D or 2D
### Author: Valay Agarawal, Dipanjali Halder, Rahul Maitra
### Date: 18 Aug, 2020

#Import required libraries. finad description of each term in the related file
import inp
import numpy as np
from inp import spin_sites, spatial_sites,dim
if dim==2:
    from inp import spin_sites_x, spin_sites_y, spatial_sites_x,spatial_sites_y
if dim==1:
    from inp import spin_sites,spatial_sites
from hubbard_hamiltonian import gen_hopping_hh, gen_onsite_hh
from multiplying_operators import multiply_two_operators, multiply_four_operators, clean_up
transformation = inp.transformation
if dim==1:
    hopping_hh=gen_hopping_hh()
if dim==2:
    hopping_hhx,hopping_hhy=gen_hopping_hh()
onsite_hh=gen_onsite_hh()

"""
To create the Hubbard Hamiltonian of the onsite term
write_onsite():
Step1 : Import onsite part of the Hubbard Hamiltonian
Step2 : Multiply four terms of the Hubbard Hamiltonian and stack up
Step3 : Club non-zero and common terms in terms of the gates


To create the Hubbard Hamiltonian of the hopping term
write_hopping():
Step1 : Check proper dimension
Step2 : Multiply two terms of  the Hubbard Hamiltonian and stack up
Step3 : Club non zero and common terms in terms of the gates
Step4 : Return 1 or 2 matrices based on dimensions
"""

#Function to write the onsite terms
def write_onsite():
    
    #Initialize the term
    t=np.zeros((0,2))
    
    for term in onsite_hh:
        
        #Multiply four operators in gate form, and return in the form of strings
        a=clean_up(multiply_four_operators(term[1],term[2],term[3],term[4]))
        
        #Remove zero coefficient terms
        a=a[a[:,0]!=0]
        
        #Multiply the  final term with the coefficient  (U) = onsite term coefficient
        a[:,0]=a[:,0]*term[0]

        #Add it to the previous terms
        t=np.vstack((t,a))
    # Add common terms and remove non zero terms
    return(clean_up(t))
def write_hopping():
    
    #First check for the dimension of the hubbard model
    if dim==1:

        #initialize the Hamiltonian with 0
        t=np.zeros((0,2))

        #Simiar process to the onsite terms, and, this time, instead of multiplying with U, multiply with t 
        for term in hopping_hh:
            a=clean_up(multiply_two_operators(term[1],term[2]))
            a[:,0]=a[:,0]*term[0]
            t=np.vstack((t,a))    
        return(clean_up(t))
    elif dim==2:
        #Similar process to above, however, we will get two matrices, one for hopping in x direction, another in Y
        t1=np.zeros((0,2))
        t2=np.zeros((0,2))
        for term in hopping_hhx:
            a=clean_up(multiply_two_operators(term[1],term[2]))
            a[:,0]=a[:,0]*term[0]
            t1=np.vstack((t1,a)) 
        for term in hopping_hhy:
            a=clean_up(multiply_two_operators(term[1],term[2]))
            a[:,0]=a[:,0]*term[0]
            t2=np.vstack((t2,a))
        return clean_up(t1), clean_up(t2)
onsite=write_onsite()
#print(onsite)
#print(write_hopping())
if dim==1:
    hop=write_hopping()
if dim==2:
    hop_x,hop_y=write_hopping()
