### Code for generating one body and two body molecular hamiltonian in terms of operators
### Author: Valay Agarawal, Dipanjali Halder, Rahul Maitra
### Date: 18 August 2020

import numpy as np
from itertools import product 
from bk_basic_functions import paritybasis, bravyikitaevbasis , inverse , pi_beta_inverse
#from bk_basic_functions import parityset , updateset , flipset , remainderset , stringchar
from inp import spin_sites

"""
Onebody_mh():
    Create one body terms in terms of the operators number and their creation and annihilation terms
Twobody_mh():
    Create two body terms in terms of the operators number and their creation and annihilation terms
    Step1 : Find all possible values of 4 body terms
    Step2 : Remove terms which are not possible due to spin violation or double excitation from the same spin orbital
    Step3 : Further remove terms and only have which will correspond to Coulomb or Exchange integrals
"""



def onebody_mh():
    one_body=[]
    for i in range (spin_sites):
        one_body.append(((i,1),(i,0)))
    return np.array(one_body)
def twobody_mh():
    a1=np.array(list(product(range(spin_sites), repeat=4)))
    two_body=[]
    for i in range (len(a1)):
        p=a1[i][0]; q=a1[i][1]; r=a1[i][2]; s=a1[i][3]
        if ((p-s)%2==0 and (q-r)%2==0 and  p!=q and r!=s):
            if (p==s and q==r and p!=q):#for direct terms
                two_body.append((0.5,(p,1),(q,1),(r,0),(s,0)))
            elif (p!=s and r!=q):#for exchange terms
                two_body.append((0.5,(p,1),(q,1),(r,0),(s,0)))
    return np.array(two_body)
#print(twobody_mh())
