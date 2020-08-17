import numpy as np
import cmath
from jw_bk import paritybasis, bravyikitaevbasis , inverse , pi_beta_inverse
from jw_bk import parityset , updateset , flipset , remainderset 
from jw_bk import write_string_of_gates_jw, write_string_of_gates_bk
from jw_bk import spin_sites


iota=complex(0,1)
X=np.array([[0,1],[1,0]])
Y=np.array([[0,-iota],[iota,0]])
Z=np.array([[1,0],[0,-1]])
I=np.array([[1,0],[0,1]]) 



def check_if_multiple(matrix1, matrix2):
    a=np.matmul(matrix1, np.linalg.inv(matrix2))
    if a[0][0]==a[1][1] and a[0][0]!=0 and a[0][1]==0:
        return a[0][0]
    else:
        return 0
   


#Take the string of gates in list form and create a hamiltonian matrix
def gen_gates_from_string(string):
    coeff=string[0]
    string=string[1]
    matrix=[]
    #count=0
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


def gen_string_from_gates(gate_matrix):
    string=[]
    sign=1
    coeff=gate_matrix[0]
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



def multiply_two_operators_jw(operator1, operator2):
    product=[]
    string1=write_string_of_gates_jw(operator1)
    string2=write_string_of_gates_jw(operator2)
    #print(string1, string2)
    for gates1 in string1:
        for gates2 in string2:     
            A1=gen_gates_from_string(gates1)
            A2=gen_gates_from_string(gates2)
            product.append(((A1[0]*A2[0]),np.matmul(A1[1],A2[1])))
    #print(product)
    gates=[]
    for gate_matrix in product:
        gates.append(gen_string_from_gates(gate_matrix))   
    return gates

def multiply_two_operators_bk(operator1, operator2):
    product=[]
    string1=write_string_of_gates_bk(operator1)
    string2=write_string_of_gates_bk(operator2)
    #print(string1, string2)
    for gates1 in string1:
        for gates2 in string2:     
            A1=gen_gates_from_string(gates1)
            A2=gen_gates_from_string(gates2)
            product.append(((A1[0]*A2[0]),np.matmul(A1[1],A2[1])))
    #print(product)
    gates=[]
    for gate_matrix in product:
        gates.append(gen_string_from_gates(gate_matrix))   
    return gates


def multiply_four_operators_jw(operator1,operator2,operator3, operator4):
    product=[]
    string1=write_string_of_gates_jw(operator1)
    string2=write_string_of_gates_jw(operator2)
    string3=write_string_of_gates_jw(operator3)
    string4=write_string_of_gates_jw(operator4)
    
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
    return np.array(gates)



def multiply_four_operators_bk(operator1,operator2,operator3, operator4):
    product=[]
    string1=write_string_of_gates_bk(operator1)
    string2=write_string_of_gates_bk(operator2)
    string3=write_string_of_gates_bk(operator3)
    string4=write_string_of_gates_bk(operator4)
    
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
    return np.array(gates)


def clean_up(gates):
    gate_unique=np.unique(gates[:,1])
    #print(len(gate_unique))
    t=np.zeros((len(gate_unique)),dtype=complex)
    for i in range(len(gate_unique)):
        for j in range(len(gates)):
            if(np.array_equal(np.array(gates[j,1]),np.array(gate_unique[i]))):
                t[i]+=gates[j,0]
    return(np.transpose(np.vstack((t,gate_unique))))
                #print(gates[j],gate_unique[i])
                



#a=write_string_of_gates_bk([2,True])
#b=write_string_of_gates_bk([2,False])
#print(a)
#print(b)
#print(multiply_two_operators_bk([2,True],[2,False]))
