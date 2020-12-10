import numpy as np
import math
# function for evaluating the pi matrices
def paritybasis(number_of_spin_sites):
    return np.triu(np.ones((number_of_spin_sites,number_of_spin_sites)))
# For multiplying the pi matrix with the occupation state          
def paritystate(x):
    parity_basis=np.matmul(paritybasis(number_of_spin_sites),x)
    parity_basis=parity_basis%2
    return(parity_basis)
# function for calculating the beta matrices
def gen_full_matrix(number_of_spin_sites):
    arr0=np.array([[1,0],[1,1]])  #arr1 stores the initial beta_2 matrix of dimension 2*2
    n=int(math.log(1<<(number_of_spin_sites-1).bit_length(),2))
    if n==1:
        return arr0
    else:
        zero1=np.zeros((2**(n-1), 2**(n-1)))
        zero2=zero1
        zero2[-1]=1
        arr=gen_full_matrix(int(2**(n-1)))
        c3=np.vstack((np.hstack((arr, zero1)), np.hstack((zero2, arr)))) 
        return c3
def brayvikitaevbasis(number_of_spin_states):
    return gen_full_matrix(number_of_spin_states)[0:number_of_spin_states,0:number_of_spin_states]
print(brayvikitaevbasis(3))
print(brayvikitaevbasis(8))
#print(gen_full_matrix(11))
def bkstate(number_of_spin_sites,x):
    bk_basis=np.matmul(brayvikitaevbasis(number_of_spin_sites),x)
    bk_basis=bk_basis%2
    return (bk_basis)
####################################################################################
# main function
#number_of_spin_sites=int(input("Enter number of spin sites: "))
number_of_spin_sites=5#int(input("Enter number of spin sites: "))
x=[]
x=np.array([1,1,0,0,1])
#for i in range (number_of_spin_sites):
#    print("Enter occupation number of site",i,end=" ")
#    x.insert(i,int(input()))
x=np.array(x)
x=x.reshape(number_of_spin_sites,1)
print("Occupation number representation:")
print(x)
print("Pi Matrix:")
pimatrix=paritybasis(number_of_spin_sites)
print(pimatrix)
print("Parity basis representation:")
print(paritystate(x))
print("Beta Matrix:")
betamatrix=brayvikitaevbasis(number_of_spin_sites)
print(betamatrix)
print("BK basis representation:")
print(bkstate(number_of_spin_sites,x))
###################################################################
'''print("Beta inverse matrix:")
betainverse=np.linalg.inv(betamatrix)
betainverse=np.array(betainverse)
betainverse=betainverse%2 # To produce the modulo 2 sum
print(betainverse)
print("The transformation matrix from BK Basis to Pi basis:")
parityset=np.matmul(pimatrix,betainverse)
parityset=parityset%2
#print(parityset)'''
