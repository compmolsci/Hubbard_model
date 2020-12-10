### Code for generating BK transformation
### Author: Dipanjali Halder, Valay Agarawal, Rahul Maitra
### Date: 18 Aug 2020



import numpy as np
import cmath
import inp
iota=complex(0,1)
X=np.array([[0,1],[1,0]])
Y=np.array([[0,-iota],[iota,0]])
Z=np.array([[1,0],[0,-1]])
I=np.array([[1,0],[0,1]]) 
spin_sites=inp.spin_sites 

"""
parity basis: Conversion matrix that converts fermionic occupation basis to qubit parity basis
bravyikitaevbasis: Conversion matrix that converts from fermionic occupation to qubit BK basis
    Step1 : Get the smallest power of 2 larger than or equal to the number of spin sites(2^{n}(>=spin_sites)
    Step2 : Recrusive 2^{n} matrix generation
    Step3 : Isolate the first spin_sites*spin_sites matrix
pi_beta_inverse: Map between parity and Bravyi Kitaev
parityset, updateset, flipset, remainderset: Generate the corresponding the sets
stringchar: Convert operator into string of gates using BK transformation
"""


#Function to define the matrix to change from fermionic basis to parity basis
def paritybasis():
    # A lower triangular matrix with all elements being one, and diagonal being zero
    return(np.tril(np.ones((spin_sites,spin_sites)),-1))
parity_basis=paritybasis()

#Funciton to define the matrix to change from fermionic basis to Bravyi-Kitaev basis
def bravyikitaevbasis():
    #First find the smallest k, such that 2^k >n (number of spin_sites)
    for k in range (0,spin_sites): 
        if(spin_sites==2**k):
            n=k+1    # n is the index number of beta matrix
        if(2**(k)<spin_sites<(2**(k+1))):
            n=k+2    # n is the index number of beta matrix
    
    #create a basic matrix of 2*2, which will initialze the bk basis
    #Form of bk basis \beta_{2^{x+1}}=[[\beta_{2^x},0],[A,\beta_{2^x}]]
    arr1=np.array([[1,0],[1,1]])  #arr1 stores the initial beta_2 matrix of dimension 2*2
    for i in range (1,n-1):         # for creating all the beta matrices 
        
        #Forming zero matrix
        zero=np.zeros((2**i,2**i))
        #Forming A matrix
        zero1=np.zeros((2**i,2**i))
        zero1[-1]=1
        
        #Putting together a total matrix
        c1=np.concatenate((arr1,zero),axis=1)
        c2=np.concatenate((zero1,arr1),axis=1)
        c3=np.concatenate((c1,c2),axis=0)
        
        #Initalizing the total matrix as base matrix for rext iteration till the desired size is reached
        arr1=c3

    #Scooping out the required terms
    return(arr1[0:spin_sites,0:spin_sites])  
bkbasis=bravyikitaevbasis()

#Function to make the inverse of bk basis
def inverse():
    return(np.linalg.inv(bkbasis))%2
inv=inverse()

def pi_beta_inverse():
    return(np.matmul(parity_basis,inv))%2
pibetainv=pi_beta_inverse()


#This gets the list of the elements of the parity set for BK basis
def parityset():
    return [np.nonzero(row)[0] for row in pibetainv]
parity=parityset()

#This gets the elements for the update set
def updateset():
    return [np.nonzero(row)[0][1:] for row in np.transpose(bkbasis)]
update=updateset()

#this gets the elements of the flip set
def flipset():
    return [np.nonzero(row)[0][:-1] for row in inv]
flip=flipset()

#This gets the elements of the remainder set
def remainderset():
    return [list(set(parity[i])-set(flip[i])) for i in range(spin_sites)]
remainder=remainderset()


#This function is to convert the operator from fermionic basis to BK basis in terms of gates
def stringchar(op):

    #The operator has two parts of information, the first part defines the target site number, and the second type gives the information about the presence of dagger
    l=op[0]
    dag=op[1]

    #We initialize the operator with two strings, both with Identity matrices
    string1=["I" for m in range(spin_sites)]
    string2=["I" for m in range(spin_sites)]

    #We now iterate over all spin sites
    for i in range(spin_sites):
        
        #If some site is in the update set, we initalize it with X gate
        for k in range (len(update[l])):
            if(i==update[l][k]):
                string1[i]="X"
                string2[i]="X"                                
        #If the site number is odd, then we add it wth parity set, else remainder set
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
        #If the operator is having a dagger, we make the strings of the form X-iY, else X+iY. The coefficients are there to since Q(\dagger) = 1/2*(X-(+)iY), so we take the coeffiient 1/2 and i outside the gate structure, and then work with only the unitary gates
        if(i==l):
            if dag:
                string1[i]="X"
                string2[i]="Y"
                a=-0.5*iota
            else:
                string1[i]="X"
                string2[i]="Y"
                a=0.5*iota
    # The output is of the form coeff1*string1, coeff2*string2
    return ((0.5,string1),(a,string2))
