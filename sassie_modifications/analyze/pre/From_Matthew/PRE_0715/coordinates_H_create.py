import numpy as np

filename = "./Python/1D3Z_mod1.pdb"
new_filename = "./Python/coordinates_H.txt"

arr = np.genfromtxt(filename, dtype='str')

arr = arr[:, 5:8].astype(np.float)

np.savetxt(new_filename, arr)
