import numpy as np
from inp import trotter_order_num,trotter_step_num
def trotter_order(A,B):
    mat=np.zeros((0,2))
    if trotter_order_num==1:
        mat=np.vstack((mat,A))
        mat=np.vstack((mat,B))
        return mat
    elif trotter_order_num==2:
        A[:,0]=A[:,0]/2
        mat=np.vstack((mat,A))
        mat=np.vstack((mat,B))
        mat=np.vstack((mat,A))
        return mat
    elif trotter_order_num==3:
        A1=A
        B1=B
        A1[:,0]=A[:,0]*(7/24) 
        mat=np.vstack((mat,A1))
        B1[:,0]=B[:,0]*(2/3)
        mat=np.vstack((mat,B1))
        A1[:,0]=A[:,0]*(3/4)
        mat=np.vstack((mat,A1))
        B1[:,0]=B[:,0]*(-2/3)
        mat=np.vstack((mat,B1))
        A1[:,0]=A[:,0]*(-1/24) 
        mat=np.vstack((mat,A1)) 
        B1[:,0]=B[:,0]*(-1) 
        mat=np.vstack((mat,B1))
        return mat
    else:
        print("Specify trotter step")

def trotter_step(A):
    mat=np.zeros((0,2))
    if isinstance(trotter_step_num,int):
        #print(A[:,0])
        A[:,0]=A[:,0]/trotter_step_num
        #print(A[:,0])
        for i in range(trotter_step_num):
            mat=np.vstack((mat,A))
        return mat
    else:
        print('Enter integer trotter step')
    
#A=np.random.random((2,2))
#B=np.random.random((2,2))
#print(A)
#print(B)
#C=trotterize(A,B)
#print(C)
