import subprocess
import sys 
import numpy as np
import pyscf
import inp
from pyscf import gto, scf
def run_pes(dist):
    process=subprocess.Popen(["./pes.sh",str(dist)],stdout=subprocess.PIPE)
    process.wait()
    a=process.communicate()[0]
    return a#.splitlines()
def get_enuc(dist):
    #Specify the molecule, according to pyscf. Input atoms, unit of distance, coordinates and basis set
    mol=gto.M(unit='ANGSTROM',atom='''H 0 0 0; H 0 0 '''+str(dist),basis='sto-3g' )
    #mol=gto.M(unit='Bohr',atom='''H 0 0 1.41; H 0 0 -1.41''',basis='sto-3g')
    nao= mol.nao_nr();nelec=mol.nelectron;occ=int(nelec/2);virt=nao-occ
    enuc=mol.energy_nuc()
    return enuc

distance=np.arange(0.5,0.7,0.05)
distance2=np.arange(0.7,0.8,0.01)
distance3=np.arange(0.8,2,0.1)

distances=np.concatenate((distance,distance2))
distances=np.concatenate((distances,distance3))

distances=np.arange(0.9,1.1,0.05)
#distances=[0.74]
dist_energy=np.zeros((len(distances),2))
#dist_energy=np.zeros((n_points,2))
i=0
#dist=0.735
for dist in distances:
#for i in range(3):
    #dist=0.5+0.1*i
    #dist=0.735
    energy=float(run_pes(dist))
    
    #print(energy)
    #energy=0
    #print(get_enuc(dist))
    dist_energy[i,0]=dist
    #print(energy
    #e1=get_enuc(dist)
    e1=inp.enuc
    print('energy after evolution', energy)
    print('nuclear_energy',e1)
    #print('bond distance', dist, 'energy', energy)
    dist_energy[i,1]=energy+e1
    print('distance and energy with nuclear terms added',dist_energy[i])
    #dist+=0.735
    i+=1
#np.save('H2_pes_8_ancillae',dist_energy)
#print(dist_energy)
