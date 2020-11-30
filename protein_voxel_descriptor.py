import numpy as np
from moleculekit.vmdviewer import getCurrentViewer, getVMDpath
from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from moleculekit.home import home
from pathlib import Path
import os


pdb_filepath = '/Users/smd4193/Desktop/2k39.pdb'

prot = Molecule(pdb_filepath)

vmd_path = getVMDpath()
curr_view = getCurrentViewer()

# print(prot)
prot = prepareProteinForAtomtyping(prot)
prot_vox, prot_centers, prot_N = getVoxelDescriptors(prot, buffer=1, voxelsize=2)
print('heho')

# x = prot_centers[:,0]
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.plot(prot_centers[:,0], prot_centers[:, 1], prot_vox[:, 2])
ax.scatter(prot_centers[:,0], prot_centers[:, 1], prot_centers[:, 2], s=np.multiply(prot_vox[:, 2], 10))
plt.show()
plt.close()

# prot.view(guessBonds=False)
viewVoxelFeatures(prot_vox, prot_centers, prot_N)
