import numpy as np
import pyscf
from math import floor
from pyscf import gto, scf


#Input the number of trotter order
trotter_order_num=1

#Input the number of trotter steps
trotter_step_num=1

#arrange='aabb'
arrange='abab'

total_time=3
#Input the time step
time_step=1
info='t_order='+str(trotter_order_num)+',t_step='+str(trotter_step_num)+',time_step='+str(time_step)+'\n'
#Input the basis of transformation
transformation = 'Jordan-Wigner';info+='JW,'
#transformation = 'Bravyi-Kitaev';info+=' JW '
print('transformation = ',transformation)

#Input the Hamiltonian you want to work in 
hamiltonian ='molecular';info+='molH '
#hamiltonian ='hubbard';info+='hubH,'


#For a Hubbard Hamiltonian
if hamiltonian == 'hubbard':
    #Specify the usage of periodic boundary condition. If true, then periodic boundary condition will be used
    boundary=False;
    if boundary==False:info+='no '
    info+='PBC,'
    
    #Specify the dimension of Hubbard Model
    dim=1

    #If the dimension is 1 (Linear Model)
    if dim==1:
        #Specify the hopping terms coefficient value
        t=1

        #Specify the onsite terms' coefficient value
        u=1

        #Number of Sites in the model
        spatial_sites=3
        spin_sites=2*spatial_sites

        #Number of electrons you want to work with
        electrons=4
        
        if spatial_sites==2:
            coupling='singlet'
            #coupling='triplet'
            initinfo=',2-'+str(coupling)
        elif spatial_sites==3:
            coupling='doublet'
            ms=0.5
            #ms=-0.5
            first_two_coupling='singlet'
            #first_two_coupling='triplet'
            initinfo=str(coupling)+',ms='+str(ms)+',first_two_coupling_'+str(first_two_coupling)
        elif spatial_sites==4:
            coupling=='singlet'
            coupling=='triplet-triplet-singlet'
        modelinfo=',sites='+str(spatial_sites)+',t/U='+str(t/u)
        print(hamiltonian,', boundary = ',boundary, ', dim =',dim,', spatial_sites =',spatial_sites,', t/U =',t/u)
    #If the dimension is 2 (2D model)
    if dim==2:

        #Specify the hopping in x direction terms' coefficient value
        tx=1
        
        #Specify the hopping in x direction terms' coefficient value
        ty=1
        
        #Specify the onsite terms' coefficient value
        u=1

        #Number of sites in x direction
        spatial_sites_x=2
        
        #Include the linearization type
        pattern='snake'

        #Number of sites in y directions
        spatial_sites_y=2

        #Number of electrons
        electrons=4
        spin_sites_x=2*spatial_sites_x
        spin_sites_y=2*spatial_sites_y
        spatial_sites=spatial_sites_x*spatial_sites_y
        spin_sites=2*spatial_sites
        print(hamiltonian, ', dim =',dim,', spatial_sites_x =',spatial_sites_x,', spatial_sites_y =',spatial_sites_y,', tx/U =',tx/u, ', ty/U = ',ty/u)
        modelinfo='sites_x='+str(spatial_sites_x)+' sites_y='+str(spatial_sites_y)+' tx/U='+str(tx/u)+' ty/U='+str(ty/u)
    info+='\n'+str(dim)+'Dim,'+modelinfo+'\n'+initinfo


if hamiltonian=='molecular':
    dist
    #Specify the molecule, according to pyscf. Input atoms, unit of distance, coordinates and basis set
    mol=gto.M(unit='ANGSTROM',atom='''H 0 0 0; H 0 0 '''+str(dist),basis='sto-3g' )
    #mol=gto.M(unit='Bohr',atom='''H 0 0 1.41; H 0 0 -1.41''',basis='sto-3g')
    nao= mol.nao_nr();nelec=mol.nelectron;occ=int(nelec/2);virt=nao-occ
    enuc=mol.energy_nuc()
    #Spatial sites will be equal to number of atomic orbitals taken
    spatial_sites=nao
    spin_sites=2*spatial_sites
    electrons=nelec
    #print('nao',nao,'occ',occ)

    #The kinetic aenergy term
    T = mol.intor('cint1e_kin_sph')
    #print(T)
    #Nuclear Repulsion term
    V = mol.intor('cint1e_nuc_sph')
    #print(V)
    #print(V)
    #Exchange terms
    v2e = mol.intor('cint2e_sph').reshape((nao,)*4)
    #print(v2e)
    #Running a scf procedure
    mf = scf.RHF(mol).run()
    #print(mol.intor('nuclear_repulsion_sph'))
            
    #extracting useful quantities
    E_hf = mf.e_tot
    hf_mo_E = mf.mo_energy
    hf_mo_coeff = mf.mo_coeff
    #print(mf.nuc)
    #Changing the fock and 2 electron repulsion matrix for to be used with other codes
    Fock = T + V
    oneelecint_mo = np.einsum('ab,ac,cd->bd',hf_mo_coeff,Fock,hf_mo_coeff) 
    oneelecint_mo[abs(oneelecint_mo) < 10**(-8)] = 0
    twoelecint_1 = np.einsum('zs,wxyz->wxys',hf_mo_coeff,v2e) 
    twoelecint_2 = np.einsum('yr,wxys->wxrs',hf_mo_coeff,twoelecint_1)
    twoelecint_3 = np.einsum('xq,wxrs->wqrs',hf_mo_coeff,twoelecint_2)
    twoelecint_mo = np.einsum('wp,wqrs->pqrs',hf_mo_coeff,twoelecint_3)
    twoelecint_mo[abs(twoelecint_mo) < 10**(-8)] = 0
    twoelecint_mo=np.swapaxes(twoelecint_mo,1,2)
    #print(twoelecint_mo)
    print(enuc)
