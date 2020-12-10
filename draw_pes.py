import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
a=np.load('H2_pes_8_ancillae.npy')
print(a[5:30])
plt.scatter(a[:,0],a[:,1])#,msize=10)
plt.xlabel('internuclear distance')
plt.ylabel('energy')
plt.xlim([0.4,2.0])
plt.title('PES $H_2$ ancilla=8, QPE,n=1,t=1')
plt.savefig('pes_h2.png')
