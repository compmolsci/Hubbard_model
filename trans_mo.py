
                   ##----------------------------------------------------------------------------------------------------------------##
                                   
                                           # Routine to calculate Hatree-Fock energy in both spinorbital and spatial #
                                                        # orbital basis and verify with pyscf routine #
                                                # Author: Soumi Tribedi, Anish Chakraborty, Rahul Maitra #
                                                                  # Date - Feb, 2019 # 

                   ##----------------------------------------------------------------------------------------------------------------##

##----------------------------------------##
       #Import important modules#
##----------------------------------------##

import numpy as np
import copy as cp
import pyscf
import inp

from pyscf import gto, scf, cc, mcscf
#from pyscf.cc import ccsd_t

##----------------------------------------------------##
       #Import important parameters from pyscf#
##----------------------------------------------------##

mol = inp.mol
#print(mol.atom)
# Obtain the number of atomic orbitals in the basis set
nao = mol.nao_nr()
# Obtain the number of electrons
nelec = mol.nelectron
#print('nelec',nelec)
# Compute nuclear repulsion energy
enuc = mol.energy_nuc()
# Compute one-electron kinetic integrals
T = mol.intor('cint1e_kin_sph')
# Compute one-electron potential integrals
V = mol.intor('cint1e_nuc_sph')
# Compute two-electron repulsion integrals (Chemists' notation)
v2e = mol.intor('cint2e_sph').reshape((nao,)*4)
#print(T.shape, V.shape, v2e.shape)
#--------------------------------------------------#
                  #Initialize#
#--------------------------------------------------#

'''if inp.iniguess == 'mcscf':
  mf1 = scf.RHF(mol).run()
  mc = mcscf.CASSCF(mf1, 4, 4)
  mc.chkfile ='/home/visitors/Anish/singlet-channel/package_spin/chk1'
  info = {'e_tot'    : mc.e_tot,
          'mo_energy': mc.mo_energy,
          'mo_coeff' : mc.mo_coeff}
  mc.kernel()
  mf = scf.RHF(mol)
  mf.init_guess = 'chk1'
  mf.kernel()'''

mf = scf.RHF(mol).run()
E_hf = mf.e_tot
hf_mo_E = mf.mo_energy
hf_mo_coeff = mf.mo_coeff
n = nelec/2

occ = int(n)
virt = nao-occ
#print(occ, virt)
#if inp.mode == 'spin':
nso = 2*nao                          # number of spin orbitals
#occ_alpha = occ                      # number of orbitals occupied by alpha electrons
#occ_beta = occ
#occ_so = occ_alpha + occ_beta
occ_so=2*occ
#virt_alpha = virt
#virt_beta = virt
#virt_so = virt_alpha + virt_beta
virt_so=2*virt

##--------------------------------------##
           #SET UP CORE FOCK#
##--------------------------------------##

Fock = T + V
print(T.shape, V.shape)
                                   ##-------------------------------------------------------------------------------------------------##
                                                                 ###     SPATIAL ORBITAL CODE   ###
                                   ##-------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------##
      #Transform the 1 electron integral to MO basis#
##--------------------------------------------------------------##

oneelecint_mo = np.einsum('ab,ac,cd->bd',hf_mo_coeff,Fock,hf_mo_coeff)
print(Fock)
##------------------------------------------------------------------------##
             #Transform 2 electron integral to MO Basis#
##------------------------------------------------------------------------##

twoelecint_1 = np.einsum('zs,wxyz->wxys',hf_mo_coeff,v2e)
twoelecint_2 = np.einsum('yr,wxys->wxrs',hf_mo_coeff,twoelecint_1)
twoelecint_3 = np.einsum('xq,wxrs->wqrs',hf_mo_coeff,twoelecint_2)
twoelecint_mo_temp = np.einsum('wp,wqrs->pqrs',hf_mo_coeff,twoelecint_3)

twoelecint_mo = np.swapaxes(twoelecint_mo_temp,1,2)  #physicist notation

##--------------------------------------------------------------##
                  #Verify integrals#
##--------------------------------------------------------------##

'''if inp.mode == 'spatial':
  E_scf_mo_1 = 0
  E_scf_mo_2 = 0
  E_scf_mo_3 = 0

  for i in range(0,n):
    E_scf_mo_1 += oneelecint_mo[i][i]

  for i in range(0,n):
    for j in range(0,n):
      E_scf_mo_2 += 2*twoelecint_mo[i][j][i][j] - twoelecint_mo[i][j][j][i]

  Escf_mo = 2*E_scf_mo_1 + E_scf_mo_2 + enuc
  print "HF energy in spatial orbital:", Escf_mo

##--------------------------------------------------------------##
                    #Create diagonal Fock matrix#
##--------------------------------------------------------------##
  
  Fock_mo = np.zeros((nao,nao))
  for i in range(0,nao):
    Fock_mo[i,i] = hf_mo_E[i]'''

                           ##----------------------------------------------------------------------------------------------------##
                                                                ###      SPIN ORBITAL CODE    ###
                           ##----------------------------------------------------------------------------------------------------##

#f inp.mode == 'spin':

##-----------------------------------------------------------------------##
        #Transform the 1 electron integral to spinorbital form#
##-----------------------------------------------------------------------##
  
'''oneelecint_SO = np.zeros((nso,nso))
oneelecint_SO[:occ,:occ] = oneelecint_mo[:occ,:occ]
oneelecint_SO[occ:nelec,occ:nelec] = oneelecint_mo[:occ,:occ]
oneelecint_SO[nelec:nelec+virt,nelec:nelec+virt] = oneelecint_mo[occ:nao,occ:nao]
oneelecint_SO[nelec+virt:,nelec+virt:] = oneelecint_mo[occ:nao,occ:nao]
print(oneelecint_mo)
print(oneelecint_SO)
##-----------------------------------------------------------------------##
        #Transform the 2 electron integral to spinorbital form#
##-----------------------------------------------------------------------##

v_so_1 = np.zeros((nso,nao,nao,nao))
v_so_2 = np.zeros((nso,nso,nao,nao))
v_so_3 = np.zeros((nso,nso,nso,nao))
v_so = np.zeros((nso,nso,nso,nso))

#---------------------------------------------------------------#
for p in range(0,nao):
    v_so_1[2*p,:,:,:] = twoelecint_mo[p,:,:,:]
    v_so_1[2*p+1,:,:,:] = twoelecint_mo[p,:,:,:]

#---------------------------------------------------------------#
for q in range(0,nao):
    v_so_2[:,2*q,:,:] = v_so_1[:,q,:,:]
    v_so_2[:,2*q+1,:,:] = v_so_1[:,q,:,:]

#---------------------------------------------------------------#
for r in range(0,nao):
    v_so_3[:,:,2*r,:] = v_so_2[:,:,r,:]
    v_so_3[:,:,2*r+1,:] = v_so_2[:,:,r,:]

#---------------------------------------------------------------#
for s in range(0,nao):
    v_so[:,:,:,2*s] = v_so_3[:,:,:,s]
    v_so[:,:,:,2*s+1] = v_so_3[:,:,:,s]
#---------------------------------------------------------------#
#  v_so_1 = None
#  v_so_2 = None
#  v_so_3 = None

##-----------------------------------------------------------------##
              #TRANSITION TO OPPOSITE SPIN COMPONENT
##-----------------------------------------------------------------##

for p in range(0,nso):                 #all alpha(even) to beta(odd) transition and vice-versa will have zero contribution
    for r in range(0,nso):
        if (p % 2) == 0 and (r%2) != 0:
            v_so[p,:,r,:] = 0.0
        if (p % 2) != 0 and (r%2) == 0:
            v_so[p,:,r,:] = 0.0
    
for q in range(0,nso):
    for s in range(0,nso):
        if (q % 2) == 0 and (s%2) != 0:
            v_so[:,q,:,s] = 0.0
        if (q % 2) != 0 and (s%2) == 0:
            v_so[:,q,:,s] = 0.0
print(twoelecint_mo)#[abs(twoelecint_mo)>10**(-10)])
#print(v_so)#[abs(v_so)>10**(-10)])
##--------------------------------------------------------------##
              #Verify integrals in spinorbital basis#
##--------------------------------------------------------------##

E_scf_mo_so_1 = 0
E_scf_mo_so_2 = 0
E_scf_mo_3 = 0
for i in range(0,nelec):
    E_scf_mo_so_1 += oneelecint_SO[i][i]

for i in range(0,nelec):
    for j in range(0,i):                                          #To avoid overcounting i>j
        E_scf_mo_so_2 += 1.00*(v_so[i][j][i][j] - v_so[i][j][j][i])

Escf_mo = E_scf_mo_so_1 + E_scf_mo_so_2 + enuc
print ("HF energy in spinorbital basis:", Escf_mo)

##--------------------------------------------------------------------------##
              #Orbital energies in spinorbital basis#
##--------------------------------------------------------------------------##
  
mo_E = np.zeros((nso))
for i in range(0,nao):
    mo_E[2*i] = hf_mo_E[i]
    mo_E[2*i+1] = hf_mo_E[i]

##--------------------------------------------------------------------------##
              #Create diagonal Fock matrix in spinorbital basis#
##--------------------------------------------------------------------------##

F_mo_SO = np.zeros((nso,nso))
for i in range(0,nso):
    F_mo_SO[i,i]=mo_E[i]


##--------------------------------------------------------------------------##
                      #Verification against pyscf#
##--------------------------------------------------------------------------##

def check_mo():
  if abs(Escf_mo - E_hf) <= 1E-6 :
    print ("MO conversion successful")
  return
print(Escf_mo,E_hf)
check_mo()

#gc.collect()'''

                          ##-----------------------------------------------------------------------------------------------------------------------##
                                                                                    #THE END#
                          ##-----------------------------------------------------------------------------------------------------------------------##

