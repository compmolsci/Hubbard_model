
#### Program to write all the Hubbard Hamiltonian terms in the terms of operators and apply them on the states
### Author: Valay Agarawal, Dipanjali Halder, Rahul Maitra
### Date: 18 Aug 2020


# Import required stuff
import numpy as np
import inp
from inp import dim
if dim==1:
    from inp import t,u,spatial_sites, spin_sites
if dim==2:
    from inp import tx,ty,u,spatial_sites_x, spin_sites_x, spatial_sites_y,spin_sites_y,pattern
    nx=spatial_sites_x
    ny=spatial_sites_y
from inp import boundary
from basic_functions import qubit_to_classical, classical_to_qubit  
from string_interconversion import apply_gates, binary_value, apply_operator, apply_string 
from multiplying_operators import multiply_two_operators, multiply_four_operators
'''
Notes:
    Calculations done assuming a UDUDUD... spin arrangement as compared to a UUU....DDD.... spin arrangement


A_{i,\sigma}^(\dagger)*A_{j,\sigma} is the form of the hopping terms
gen_hopping_hh()
Generates all the terms of hopping term of Hubbard Hamiltonian in the form of operator identifiers
Possibilites: 
Dimension = 1 or 2
Periodic Boundary Conditions: True or False
Step1 : Check for correct Boundary Condition and Dimension
Step2 : Write the Hamiltonian in the form of only spin sites, where 1st term will be creation, second term will be annihilation operator
Step3 : Change the spin site only hamiltonian to include the information of the creation/annihilation
gen_onsite_hh(operators):
Generates all the terms of onsite term of Hubbard Hamiltonian in the form of operator identifiers
A_{i,\sigma_1}^(\dagger)*A_{i,\sigma_1}*A_{i\sigma_2}^{\dagger}*A_{i\sigma_2}
Here, the dimension does not matter, so repeat Step 2 and Step 3 of gen_hopping_hh()

apply_onsite_hh(onsite_hh,all_states, all_states_binary):
Applies all the onsite terms of the Hubbard Hamiltonian and gives out a coefficient matrix

apply_hopping_hh(onsite_hh,all_states, all_states_binary):
Applies all the hopping terms of the Hubbard Hamiltonian and gives out a coefficient matrix

'''




spatial_sites=inp.spatial_sites
spin_sites=2*spatial_sites
#Generate hopping hamiltonian part in terms of operators
def gen_hopping_hh():
    #First write operator indices, and then add the information of creation/annihilation and the coefficients to be multiplied
    
    #check for appropriate dimension
    if dim==1:
        hopping_hh=[]
        #check for appropriate boundary condition
        if boundary==True:
            
            #Include A->B and B->A hopping, where index of A < Index of B. And then iterate over all values of A
            for i in range(spin_sites):
                hopping_hh.append(([i,(i+2)%spin_sites]))
                hopping_hh.append(([(i+2)%spin_sites,i]))
            hopping_hh=np.unique(np.array(hopping_hh),axis=0)
        else:
            #Do the same as above, except when A->B hopping would mean going out of the range of spin sites
            for i in range(spin_sites-2):
                hopping_hh.append(([i,i+2]))
                hopping_hh.append(([i+2,i]))
        hopping=[]
        for i in range(len(hopping_hh)):

            #Include the information of hopping coefficients and dagger
            hopping.append((-t,[hopping_hh[i][0],1],[hopping_hh[i][1],0]))
        return np.array(hopping)
    
    
    elif dim==2:
        #For snake pattern, we follow x major form, 
        if pattern=='snake':
            hopping_hhx=[]
            hopping_hhy=[]
            if boundary==True:
                for i in range(nx):
                    for j in range(ny):  
                        for k in range(2):
                            ''' The terms are of the form A->B, B->A, A->C, C->A. 
                            A_x = Index in x direction
                            A_y = Index in y direction
                            B_x=A_x+1
                            B_y=A_y+1
                            C_x=A_x
                            C_y=C_y+1
                            '''
                            hopping_hhx.append(([2*i*ny+2*j+k,2*ny*((i+1)%nx)+2*j+k]))
                            hopping_hhx.append(([2*((i+1)%nx)*ny+2*j+k,2*ny*i+2*j+k]))
                            hopping_hhy.append(([2*i*ny+2*j+k,2*ny*i+2*((j+1)%ny)+k]))
                            hopping_hhy.append(([2*i*ny+2*((j+1)%ny)+k,2*ny*i+2*j+k]))
                hopping_hhx=np.unique(np.array(hopping_hhx),axis=0)
                hopping_hhy=np.unique(np.array(hopping_hhy),axis=0)
            else:
                if pattern=='snake':
                    for i in range(nx):
                        for j in range(ny):
                            for k in range(2):
                                ''' The terms are of the form A->B, B->A, A->C, C->A.
                                    A_x = Index in x direction
                                    A_y = Index in y direction
                                    B_x=A_x+1
                                    B_y=A_y+1
                                    C_x=A_x
                                    C_y=C_y+1
                                '''
                                #Except for the last term, rest is the same as without boundary conditions. 
                                if i!=nx-1:
                                    hopping_hhx.append(([2*i*ny+2*j+k,2*ny*(i+1)+2*j+k]))
                                    hopping_hhx.append(([2*ny*(i+1)+2*j+k,2*ny*i+2*j+k]))
                                if j!=ny-1:
                                    hopping_hhy.append(([2*i*ny+2*j+k,2*ny*i+2*(j+1)+k]))
                                    hopping_hhy.append(([2*i*ny+2*(j+1)+k,2*ny*i+2*j+k]))
        hoppingx=[]
        hoppingy=[]
        # Add the information of hopping coefficient and creation/annihilation to both hopping_hhx and hopping_hhy 
        for i in range(len(hopping_hhy)):
            hoppingy.append((-tx,[hopping_hhy[i][0],1],[hopping_hhy[i][1],0]))
        for i in range(len(hopping_hhx)):
            hoppingx.append((-ty,[hopping_hhx[i][0],1],[hopping_hhx[i][1],0]))
        return np.array(hoppingx),np.array(hoppingy)



#Generate onsite part of Hubbard Hamitonian in terms of operators
def gen_onsite_hh():
    onsite=[]
    #check for appropriate dimension
    if dim==1:
        for i in range(spatial_sites):
            onsite.append((-u,[2*i,1],[2*i,0],[2*i+1,1],[2*i+1,0]))
    elif dim==2:
        if pattern=='snake':
            for i in range(nx):
                for j in range(ny):
                    onsite.append((-u,[2*nx*i+2*j,1],[2*nx*i+2*j,0],[2*nx*i+2*j+1,1],[2*nx*i+2*j+1,0]))
    return np.array(onsite)
#print(gen_onsite_hh())

#Apply onsite part of hubbard hamiltonian
def apply_onsite_hh(onsite_hh,all_states, all_states_binary):
    matrix=np.zeros((len(all_states),len(all_states)),dtype=complex)
    #print(matrix)
    for x in range(len(onsite_hh)):
        for inp_state in all_states:
            coeff=inp_state[0]
            inp_state=inp_state[1]
            #print('inp_state',inp_state)
            string_matrix=multiply_four_operators(onsite_hh[x,0],onsite_hh[x,1],onsite_hh[x,2],onsite_hh[x,3])
            #print(operator_matrix.size)
            for i in range(len(string_matrix)):
                final_coeff, final_state = apply_string(coeff, inp_state, string_matrix[i])
                if final_state!=0:
                    t1=np.where(all_states_binary == binary_value(final_coeff,final_state))[0][0]
                    t2=np.where(all_states_binary == binary_value(coeff, inp_state))[0][0]
                    matrix[t1,t2] +=final_coeff
    return matrix

#Apply hopping part of hubbard hamiltonian
def apply_hopping_hh(hopping_hh, all_states, all_states_binary):
    matrix=np.zeros((len(all_states),len(all_states)), dtype=complex)
    for x in range(len(hopping_hh)):
        for inp_state in all_states:
            coeff=inp_state[0]
            inp_state=inp_state[1]
            string_matrix=multiply_two_operators(hopping_hh[x,0],hopping_hh[x,1])
            for i in range(len(string_matrix)):
                final_coeff,final_state=apply_string(coeff,inp_state,string_matrix[i])
                if final_state!=0:
                    t1=np.where(all_states_binary == binary_value(final_coeff,final_state))
                    t2=np.where(all_states_binary == binary_value(coeff, inp_state))
                    matrix[t1,t2]+=final_coeff
    return matrix
#print(gen_hopping_hh())
