import numpy as np
import cmath
from itertools import product
from jw_bk import spin_sites
from Products import multiply_two_operators_jw, multiply_two_operators_bk, multiply_four_operators_jw, multiply_four_operators_bk
from Products import clean_up

#Hopping terms
print("Hopping terms:")
a=np.array(list(product(range(spin_sites), repeat=2)))
hopping=[]
for i in range (a.shape[0]): #Generates all the combinations of hopping term (doesn't consider periodic boundary condition)
    if ((a[i][0]-a[i][1]==2) or (a[i][0]-a[i][1]==-2)):
        hopping.append(a[i])

        
if (spin_sites!=4):  # For including the periodic boundary conditions
    hopping.append([spin_sites-1,1])
    hopping.append([spin_sites-2,0])
    hopping.append([1,spin_sites-1])
    hopping.append([0,spin_sites-2])
    
hopping=np.array(hopping)  
print(hopping) #prints all possible combinations of hopping terms

m=[]  # m stores the multiplied string for hopping part of the Hamiltonian
for i in range (hopping.shape[0]):
    p=hopping[i][0]
    q=hopping[i][1]
    m.append(multiply_two_operators_jw([p,True],[q,False]))
print(m)



#On-site terms
print("Onsite terms:")
b=np.array(list(product(range(spin_sites),repeat=4)))
onsite=[]
for i in range (b.shape[0]):
    if(b[i][0]==b[i][1] and b[i][2]==b[i][3] and b[i][2]-b[i][1]==1 and b[i][0]%2==0):
        onsite.append(b[i])
onsite=np.array(onsite)
print(onsite)

m1=[]  # m1 stores the multiplied string for onsite part of the Hamiltonian
for i in range (onsite.shape[0]):
    p=onsite[i][0]
    q=onsite[i][1]
    r=onsite[i][2]
    s=onsite[i][3]
    m1.append(multiply_four_operators_jw([p,True],[q,False],[r,True],[s,False]))
print(m1)





    
    
    