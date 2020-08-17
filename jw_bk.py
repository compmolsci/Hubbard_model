import numpy as np
import cmath
iota=complex(0,1)
X=np.array([[0,1],[1,0]])
Y=np.array([[0,-iota],[iota,0]])
Z=np.array([[1,0],[0,-1]])
I=np.array([[1,0],[0,1]]) 
spin_sites= 4 # Make it 4 to study hydrogen molecular Hamiltonian/two site Hubbard model
def paritybasis():
    return(np.tril(np.ones((spin_sites,spin_sites)),-1))
    
# function for calculating the beta matrices
def bravyikitaevbasis():
    for k in range (0,spin_sites): 
        if(spin_sites==2**k):
            n=k+1    # n is the index number of beta matrix
        if(2**(k)<spin_sites<(2**(k+1))):
            n=k+2    # n is the index number of beta matrix
      
    arr1=np.array([[1,0],[1,1]])  #arr1 stores the initial beta_2 matrix of dimension 2*2
    for i in range (1,n-1):         # for creating all the beta matrices 
        zero=np.zeros((2**i,2**i))
        zero1=np.zeros((2**i,2**i))
        zero1[-1]=1
        c1=np.concatenate((arr1,zero),axis=1)
        c2=np.concatenate((zero1,arr1),axis=1)
        c3=np.concatenate((c1,c2),axis=0)
        arr1=c3
    return(arr1[0:spin_sites,0:spin_sites])  
       
# Evaluating beta inverse
def inverse():
    return(np.linalg.inv(bravyikitaevbasis()))%2

# Evaluating pi*beta inverse
def pi_beta_inverse():
    return(np.matmul(paritybasis(),inverse()))%2

# Parity set 
def parityset():
    return [np.nonzero(row)[0] for row in pi_beta_inverse()]

# Update set
def updateset():
   #print('update_set')
    #print(np.transpose(bravyikitaevbasis()))
    return [np.nonzero(row)[0][1:] for row in np.transpose(bravyikitaevbasis())]

# Flip set
def flipset():
    #print('flip_set')
    #print(inverse())
    return [np.nonzero(row)[0][:-1] for row in inverse()]

# Remainder set
def remainderset():
    parity=parityset()
    flip=flipset()
    return [list(set(parity[i])-set(flip[i])) for i in range(spin_sites)]



# Function for printing gate sequence of creation/annihilation operator (character form)
def write_string_of_gates_bk(op_def):
    l=op_def[0]
    op=op_def[1]
    string1=["I" for m in range(spin_sites)]
    string2=["I" for m in range(spin_sites)]
    update=updateset()
    parity=parityset()
    remainder=remainderset()
    for i in range(spin_sites):
        for k in range (len(update[l])):
            if(i==update[l][k]):
                string1[i]="X"
                string2[i]="X"                                
        if(l%2==0):
            for k in range (len(parity[l])):
                if(i==parity[l][k]):
                    string1[i]="Z"
                    string2[i]="Z"
        else:
            for k in range (len(parity[l])):
                if(i==parity[l][k]):
                    string1[i]="Z"
            for k in range (len(remainder[l])):
                if(i==remainder[l][k]):
                    string2[i]="Z"            
        if(i==l):
            if (op==True): #True means dagger/creation
                string1[i]="X"
                string2[i]="Y"
                a=-0.5*iota
            else: #False means annihilation
                string1[i]="X"
                string2[i]="Y"
                a=0.5*iota
    return ((0.5,string1),(a,string2))

def write_string_of_gates_jw(op_def):
    l=op_def[0]
    op=op_def[1]
    string1=[]
    string2=[]
    for k in range(spin_sites):
        if (k<l):
            string1.append('Z')
            string2.append('Z')
        elif (k==l):
            if op==True:#True means dagger/creation operator
                string1.append('X')
                string2.append('Y')
                a=-0.5*iota
            else: #False means annihilation
                string1.append('X')
                string2.append('Y')
                a=0.5*iota
        else:
            string1.append('I')
            string2.append('I')
    return((0.5,string1),(a,string2))

#print("Jordan-Wigner:")
#print(write_string_of_gates_jw([2,False]))
#print("Bravyi-Kitaev:")
#print(write_string_of_gates_bk([2,False]))
    


